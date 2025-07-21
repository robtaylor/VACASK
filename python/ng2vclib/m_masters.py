
import re
from .generators import traverse
from .exc import ConverterError

class MastersMixin:
    pat_paramassign = re.compile(r'^[A-Za-z_][A-Za-z0-9_]*=')

    def split_instance(self, line, in_sub=None):
        """
        Splits an instance into fragments. 

        Returns 
        * name
        * a list of parts 
        * index of first model in the list of fragments, None if no model found. 

        Tracks model usage. 
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
                    key = (part, in_sub)
                    if key not in self.data["model_usage"]:
                        self.data["model_usage"][key] = set()
                    self.data["model_usage"][key].add(in_sub)
                    break
            # Global models
            if None in self.data["models"] and part in self.data["models"][None]:
                # Found model
                mod_index = ndx
                key = (part, None)
                if key not in self.data["model_usage"]:
                    self.data["model_usage"][key] = set()
                self.data["model_usage"][key].add(in_sub)
                break
        
        return name, parts, mod_index
    
    def preprocess_instance(self, lnum, lws, line, eol, annot, in_sub):
        """
        First stage of instance procesing. 
        
        Split card into name and words list, find model. 
        Store name, list of words, and model index in annotations. 
        
        Track used models. 
        """
        name, parts, mod_index = self.split_instance(line, in_sub)
        annot["name"] = name
        annot["words"] = parts
        annot["mod_index"] = mod_index
        annot["mod_name"] = parts[mod_index] if mod_index is not None else None

    def collect_masters(self):
        """
        Colects defined subcircuits and models from deck.  
        """
        deck = self.data["deck"]

        in_sub = None

        # Pass 1 - collect models
        # Pass 2 - preprocess instances, construct model usage table
        for passno in  [1, 2]:
            for history, line, depth, in_control_block in traverse(deck, depth=self.cfg.get("process_depth", None)):
                if in_control_block and passno==1:
                    # In control block, pass 1, look for pre_osdi
                    lnum, lws, l, eolc, annot = line
                    l1 = l.strip()
                    if l1.startswith("*"):
                        l1 = l1[1:].strip()
                        if l1.startswith("pre_osdi"):
                            # Have a pre_osdi, add to list of files
                            osdi_file = l1[8:].strip()
                            self.data["osdi_loads"].add(osdi_file)
                    continue

                lnum, lws, l, eolc, annot = line

                if isinstance(eolc, tuple):
                    # If eolc is a tuple this is an .include/.lib line
                    # We do not process those. 
                    continue
                elif l.startswith("*"):
                    # Comment
                    continue
                elif len(l)==0:
                    # Empty line
                    continue
                elif l.startswith(".subckt"):
                    parts = l.split(" ")
                    in_sub = parts[1]

                    if passno==1:
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
                elif l.startswith(".ends"):
                    # End of subcircuit
                    in_sub = None
                elif l.startswith(".model"):
                    if passno==1:
                        # Detect model
                        parts = l.split(" ")
                        name = parts[1]
                        mtype = parts[2]
                        params = parts[3:]
                        
                        annot["name"] = name

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
                elif not l.startswith("."):
                    if passno==2:
                        # Not a dot command, must be an instance
                        # Preprocess it
                        try:
                            self.preprocess_instance(*line, in_sub)
                        except ConverterError as e:
                            raise ConverterError(str(e), history, lnum)

            



                
