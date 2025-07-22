import re

from .generators import traverse, format_history
from .exc import ConverterError

pat_siprefix = re.compile(r'\b(\d+\.\d*|\.\d+|\d+\.|\d+)(meg|g|t|mil)\b')

prefix_map = {
    "meg": "M", 
    "g": "G", 
    "t": "T"
}

def si_replace_worker(match):
    """
    Converts SI prefixes based on a rexexp match. 
    xmeg -> xM
    xg -> xG
    xt -> xT
    xmil -> (x*25.4e-6)
    """
    num = match.group(1)
    prefix = match.group(2)
    if prefix=="mil":
        return "("+num+"*25.4e-6)"
    else:
        return num+prefix_map[prefix]

def si_replace(expr):
    """
    Converts SI prefixes to VACASK syntax. 
    """
    return pat_siprefix.sub(si_replace_worker, expr)

class OutputMixin:
    def format_extra_params(self, extras, lws, indent):
        """
        Format extra parameters.

        For now, we do not care about line length. 
        """
        intxt = lws + (" "*indent)

        txt = ""
        first = True
        for name, value in extras.items():
            if not first:
                txt += " "
            first = False
            txt += name+"="+value
        
        return intxt+txt
    
    def format_value(self, valuestr):
        """
        Formats a value. Removes curly braces. Changes SI prefixes. 
        """
        if valuestr[0]=="{":
            valuestr = valuestr[1:-1]
        
        # Handle SI prefixes
        return si_replace(valuestr)

    def format_params(self, params, atcol=0, indent=2):
        """
        Format parameters, indent them, make sure lines are not too long. 

        *lws* is the leading whitespace string. 

        *atcol* is the column number where we start to dump parameters. 

        Does not indent. Uses *indent* for computing when to split a line. 

        Returns a string, a flag indicating string has multiple lines, 
        and a flag indicating that the formatted output has more then one 
        line. 
        """
        # Strip braces from parameter values, replace SI units
        # Compute full size
        newpar = []
        full_size = 0
        for name, value in params:
            value = self.format_value(value)
            newpar.append((name, value))

            full_size += len(name)+1+len(value)+1
        params = newpar

        # Do we need to split the line 
        need_split = atcol+full_size>self.cfg["columns"]

        split = False
        line = ""
        fragment = ""
        first = True
        for name, value in params:
            # Remove curly braces
            if value[0]=="{":
                value = value[1:-1]
            
            # Handle SI prefixes
            value = si_replace(value)

            # Do we need to start a new line
            if indent+len(fragment)+len(name)+1+len(value)>self.cfg["columns"]:
                # New line
                if not first:
                    line += "\n"
                    split = True
                else:
                    first = False
                line += fragment
                fragment = ""
            
            # Space as separator
            if len(fragment)>0:
                fragment += " "
            
            # Add parameter
            fragment += name + "=" + value
        
        # Flush last fragment
        if len(fragment)>0:
            if not first:
                line += "\n"
                split = True
            line += fragment
        
        return line, need_split, split

    def indent(self, txt, indent):
        """
        Indents text by *indent* spaces. 
        """
        istr = " "*indent
        return "\n".join([istr+l for l in txt.split("\n")])

    def format_subckt_params(self, params, lws):
        """
        Format subcircuit parameters. 

        One parameter per line. 
        """
        txt = ""
        first = True
        for name, value in params:
            value = si_replace(value)
            if not first:
                txt += "\n"
            first = False
            txt += lws + "parameters "+name+"="+value

        return txt

    def split_params(self, params, handle_m=False):
        """
        Splits a list of parameter assignments into a list of 
        (names, value) pairs. 
        """
        psplit = []
        for p in params:
            split = p.split("=")
            if len(split)!=2:
                raise ConverterError("Malformed parameter '"+p+"'.")
            
            # Handle mfacto
            if handle_m:
                if p[0]=="m" or p[0]=="_mfactor":
                    split = ( "$mfactor", p[1] )
            
            psplit.append(split)
        
        return psplit

    def remove_params(self, params, to_remove=set()):
        """
        Removes all parameters listed in set *to_remove*. 
        Can handle unsplit and split parameters. 
        """
        out = []
        for part in params:
            if isinstance(part, list):
                if part[0] in to_remove:
                    continue
            elif part in to_remove:
                    continue
            out.append(part)
        return out
    def merge_vectors(self, params, vecnames=set()):
        """
        Merges vector parameters into one string. 

        Assumes *params* is a list of unsplit parameters. 
        """
        out = []
        for ndx in range(len(params)):
            part = params[ndx]
            sp = part.split("=", 1)
            if len(sp)>1 and sp[0] in vecnames:
                # Start merging
                merged = ""
                while ndx<len(params): 
                    merged += params[ndx]
                    if params[ndx].strip()[-1]!=",":
                        break
                out.append(merged)
            else:
                out.append(part)
        
        return out

    def vacask_file(self):
        """
        Returns VACASK file as a list of lines. 
        """
        deck = self.data["deck"]
        out = []
        first = True
        in_sub = None
        target_depth = self.cfg.get("output_depth", None)
        for history, line, depth, in_control_block in traverse(deck, depth=target_depth):
            # Skip control block
            if in_control_block:
                continue

            lnum, lws, l, eolc, annot = line

            # Special handling for title
            if first:
                first = False
                if self.data["is_toplevel"]:
                    out.append(self.data["title"])
                    continue
            
            # Handle various lines
            if len(l)==0:
                # Empty line
                out.append(lws+l)
            elif l.startswith("*"):
                # Comment
                out.append(lws+"//"+l)
            elif l.startswith(".model"):
                # Model, output if it is used or output is forced
                if (
                    self.cfg.get("all_models", False) or 
                    len(self.data["model_usage"].get((annot["name"], in_sub), set()))>0
                ):
                    out.append(self.process_model(lws, l, eolc, in_sub))
            elif l.startswith(".include"):
                # Include
                if target_depth is None or depth<target_depth:
                    # Not at the deepest level yet, 
                    # output a comment containing original .include
                    out.append(lws+"// "+l)
                else:
                    name, section, subdeck = eolc
                    out.append(lws+"include \""+name+"\"")
            elif l.startswith(".lib"):
                # Lib
                name, section, subdeck = eolc
                if name is None:
                    # Section start marker
                    out.append(lws+"lib "+section)
                else:
                    # Library section include
                    if target_depth is None or depth<target_depth:
                        # Not at the deepest level yet, 
                        # output a comment containing original .lib
                        out.append(lws+"// "+l)
                    else:
                        name, section, subdeck = eolc
                        out.append(lws+"include \""+name+"\" "+section)
            elif l.startswith(".endl"):
                # End of section marker
                out.append(lws+"endl")
            elif l.startswith(".subckt"):
                name = l.split(" ")[1]
                in_sub = name
                terminals, params = self.data["subckts"][name]
                vcline = lws+"subckt "+name+"("
                vcline += " ".join(terminals)
                vcline +=")"
                out.append(vcline)
                if len(params)>0:
                    out.append(self.format_subckt_params(params, lws))
            elif l.startswith(".ends"):
                in_sub = None
                out.append(lws+"ends")
            elif l.startswith(".param"):
                txt = l.replace(".param", "parameters")
                out.append(lws+txt)
            elif l.startswith(".if"):
                txt = l.replace(".if", "@if")
                out.append(lws+txt)
            elif l.startswith(".elseif"):
                txt = l.replace(".elseif", "@elseif")
                out.append(lws+txt)
            elif l.startswith(".else"):
                txt = l.replace(".else", "@else")
                out.append(lws+txt)
            elif l.startswith(".endif"):
                txt = l.replace(".endif", "@end")
                out.append(lws+txt)
            elif l.startswith(".end"):
                # Done. 
                break
            else:
                # Instance
                if l[0]>="a" and l[0]<="z":
                    method=getattr(self, "process_instance_"+l[0], None)
                    if method is None:
                        raise ConverterError("Dont' know how to process instances of type '"+l[0]+"'.", history, lnum)
                    try:
                        txt = method(lws, l, eolc, annot, in_sub)
                    except ConverterError as e:
                        raise ConverterError(str(e), history, lnum)
                    out.append(txt)
                
        # Dump load statements and default models. 
        # This will happen only in the toplevel
        m = self.load_statements()
        if len(m)>0:
            out.append("")
            out.append("// Modules")
            out.extend(m)
        
        m = self.default_models()
        if len(m)>0:
            out.append("")
            out.append("// Default models")
            out.extend(m)
        
        return out
                



            

            

