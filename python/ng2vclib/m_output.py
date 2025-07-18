import re

from .generators import traverse, format_history

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
        for extra in extras:
            if not first:
                txt += " "
            first = False
            name, value = extra
            txt += name+"="+value
        
        return intxt+txt
    
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
            # Remove curly braces
            if value[0]=="{":
                value = value[1:-1]
            
            # Handle SI prefixes
            value = si_replace(value)

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

    def split_instance(self, line, in_sub=None):
        """
        Splits an instance into fragments. 

        Returns 
        * name
        * a list of parts 
        * index of first model in the list of fragments, None if no model found
        """
        parts = line.split()
        name = parts[0]
        parts = parts[1:]

        # Find first model
        mod_index = None
        for ndx, part in enumerate(parts):
            # Local models
            if in_sub is not None:
                if in_sub in self.data["models"] and part in self.data["models"][in_sub]:
                    # Found model
                    mod_index = ndx
                    break
            # Global models
            if None in self.data["models"] and part in self.data["models"][None]:
                # Found model
                mod_index = ndx
                break
        
        return name, parts, mod_index
    
    def vacask_file(self):
        """
        Returns VACASK file as list of lines. 
        """
        deck = self.data["deck"]
        out = []
        first = True
        in_sub = None
        for history, line, in_control_block in traverse(deck, recursive=self.cfg.get("recursive_process", False)):
            # Skip control block
            if in_control_block:
                continue

            lnum, lws, l, eolc = line

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
                # Model
                out.append(self.process_model(lws, l, eolc, in_sub))
            elif l.startswith(".include"):
                if self.cfg["flat"]:
                    out.append(lws+"// "+l)
                    # TODO
                    pass
                else:
                    subfile, name = eolc
                    out.append(lws+"include \""+name+"\"")
            elif l.startswith(".lib"):
                if self.cfg["flat"]:
                    out.append(lws+"// "+l)
                    # TODO
                    pass
                else:
                    subfile, name, section = eolc
                    out.append(lws+"include \""+name+"\" "+section)
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
            else:
                # Instance
                if l[0]>="a" and l[0]<="z":
                    method=getattr(self, "process_instance_"+l[0], None)
                    if method is None:
                        raise Exception(
                            format_history(history, lnum)+
                            "\n  Dont' know how to process instances of type '"+l[0]+"'."
                        )
                    out.append(method(lws, l, eolc, in_sub))
                
            
        return out
                



            

            

