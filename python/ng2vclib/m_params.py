import re

from .exc import ConverterError

pat_siprefix = re.compile(r'\b(\d+\.\d*|\.\d+|\d+\.|\d+)(meg|g|t|mil)\b')

pat_temper = re.compile(r"\btemper\b")

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


class ParamsMixin:
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

        Treats boolean parameters as <param>=1. 

        If requested, renames m parameter. 
        """
        psplit = []
        for p in params:
            split = p.split("=")
            if len(split)==1:
                # Handle booleans
                split = (split[0], "1")
            elif len(split)!=2:
                raise ConverterError("Malformed parameter '"+p+"'.")
            
            # Handle mfactor
            if handle_m:
                if split[0]=="m" or split[0]=="_mfactor":
                    split = ( "$mfactor", split[1] )
            
            psplit.append(split)
        
        return psplit

    def remove_params(self, params, to_remove=set()):
        """
        Removes all parameters listed in set *to_remove*. 
        Can handle unsplit and split parameters. 
        """
        out = []
        for part in params:
            if isinstance(part, list) or isinstance(part, tuple):
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

    def process_expressions(self, params):
        """
        Processes expressions, replace 
          temper -> $temp
        """
        pout = []
        for p, e in params:
            e = re.sub(pat_temper, "$temp", e)
            pout.append((p, e))
        
        return pout

    def process_instance_params(self, params, insttype, handle_m=False):
        """
        Processes instance parameters. 

        Merges vectors, removes unneeded parameters. 

        If requested, renames m parameter. 

        Returns split parameters. 
        """
        # Merge vector parameters
        vecnames = self.cfg["merge_vector_instance_params"].get(insttype, None)
        if vecnames is not None:
            params = self.merge_vectors(params, vecnames)

        # Split parameters, default boolean parameters, handle m
        psplit = self.split_params(params, handle_m=handle_m)

        # Remove unneeded parameters
        vecnames = self.cfg["remove_instance_params"].get(insttype, None)
        if vecnames is not None:
            psplit = self.remove_params(psplit, vecnames)
        
        # Process expressions
        psplit = self.process_expressions(psplit)

        return psplit
    
    def process_terminals(self, terminals):
        """
        Processes terminals. For now replaces ! with _. 
        """
        return [ t.replace("!", "_") for t in terminals]
