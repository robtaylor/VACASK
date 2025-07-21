from .exc import ConverterError

class InstancePassiveMixin:
    def process_instance_r(self, lws, line, eol, annot, in_sub):
        """
        Process R instance (resistor). 

        rname p n <r=value>|<value> [<p1=value1> <p2=value2> ...]
        rname p n <model> [<p1=value1> <p2=value2> ...]
        rname p n <value> <model> [<p1=value1> <p2=value2> ...]
        """
        name = annot["name"]
        parts = annot["words"]
        mod_index = annot["mod_index"]
        model = annot["mod_name"]
        
        terminals = parts[:2]

        if model is None:
            # No model specified
            model = self.cfg["default_model_prefix"]+"r"
            self.data["default_models_needed"].add("r")

            # Check if part 3 is not a parameter assignment
            if "=" not in parts[2]:
                # Part 3 is the resistance, the rest are parameter assignments
                params = parts[3:]

                # Add value as first parameter assignment
                psplit = [("r", self.format_value(parts[2]))]
            else:
                # Part 3 is a parameter assignment
                params = parts[2:]

                # No parameter assignments yet
                psplit = []
            
        else:
            # Have model
            if mod_index==2:
                # Second entry, immediately after terminals, no value
                psplit = []
                model = parts[mod_index]
                params = parts[3:]
            elif mod_index==3:
                # Model is 4th entry, 3rd entry must be a value
                psplit = [("r", self.format_value(parts[2]))]
                model = parts[mod_index]
                params = parts[4:]
            else:
                # Don't know how to handle
                raise ConverterError("Cannot handle model at position "+str(mod_index+1)+".")

        # Add parameter assignments
        psplit = psplit + split_params(params)
                
        txt = lws + name + " (" + (" ".join(terminals))+") "+model+" "

        if len(psplit)>0:
            fmted, need_split, split = self.format_params(psplit, len(txt))
            if need_split or split:
                fmted = self.indent(fmted, len(lws)+2)
                txt += "(" + eol + "\n" + fmted
                txt += "\n" + lws + ")"
            else:
                txt += " " + fmted
        
        return txt
