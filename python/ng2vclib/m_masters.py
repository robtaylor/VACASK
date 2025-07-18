import re
from .generators import traverse

class MastersMixin:
    pat_paramassign = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*=')

    def collect_masters(self):
        """
        Colects defined subcircuits and models from deck.  
        """
        deck = self.data["deck"]

        in_sub = None

        for history, line, in_control_block in traverse(deck, recursive=self.cfg.get("recursive_process", False)):
            # Skip control block
            if in_control_block:
                continue

            lnum, lws, l, eolc = line

            # If eolc is a tuple this is an .include/.lib line
            # We do not process those. 
            if isinstance(eolc, tuple):
                continue

            # Detect subckt, nested subck definitions are not allowed in Ngspice
            if l.startswith(".subckt"):
                parts = l.split(" ")
                in_sub = parts[1]
                
                # Find first parameter
                pstart = 2
                for pstart, s in enumerate(parts):
                    if self.pat_paramassign.match(s):
                        break
                
                # Get terminals
                terminals = parts[2:pstart]

                # Get parameters, split into name and value
                sub_params = [ p.split("=") for p in parts[pstart:]]

                # Store
                self.data["subckts"][in_sub] = (terminals, sub_params)
                
            # End of subcircuit
            if l.startswith(".ends"):
                in_sub = None

            # Detect model
            if l.startswith(".model"):
                parts = l.split(" ")
                name = parts[1]
                mtype = parts[2]
                params = parts[3:]

                # Split params into name-value pairs
                params = [p.split("=") for p in params]

                # Get information on model type
                level = None
                version = None
                builtin = False
                if mtype in self.cfg["type_map"]:
                    # Builtins needs special handling
                    builtin = True
                    extra_params, family, remove_level, remove_version = self.cfg["type_map"][mtype]
                    
                    # Collect level and version and remove them if requested
                    pnew = []
                    for pname, pval in params:
                        if pname=="level":
                            level = int(pval)
                            if remove_level:
                                continue
                        elif pname=="version":
                            version = int(pval)
                            if remove_version:
                                continue
                        else:
                            pnew.append((pname, pval))
                    
                    params = pnew

                # Add to list of models
                if in_sub not in self.data["models"]:
                    self.data["models"][in_sub] = {}
                self.data["models"][in_sub][name] = (
                    builtin, mtype, level, version, params
                )


            
