
class InstancePassiveMixin:
    def process_instance_r(self, lws, line, eol, in_sub):
        """
        Process R instance (resistor).
        """
        name, parts, mod_index = self.split_instance(line, in_sub)
        
        terminals = parts[:2]

        if mod_index is None:
            # No model specified
            model = "default_r"
            # Check if part 3 is not a paremeter assignment
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
                raise Exception("Cannot handle model at position "+str(mod_index+1)+".")

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
