import re

class MastersMixin:
    pat_paramassign = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*=')

    def collect_masters(self, lines, depth=0):
        """
        Colects defined subcircuits and models. 
        """
        in_sub = None

        for ll in lines:
            lws, l, eol = ll
            # Recursively process lines in included files
            if isinstance(eol, list):
                self.collect_masters(eol[0], depth+1)
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
                self.data["subckts"][in_sub] = (depth, terminals, sub_params)
                
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
                    depth, builtin, mtype, level, version, params
                )


            
