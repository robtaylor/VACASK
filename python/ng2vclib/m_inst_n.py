
class InstanceNMixin:
    def process_instance_n(self, lws, line, eol, in_sub):
        """
        Process N instance (OSDI device).
        """
        name, parts, mod_index = self.split_instance(line, in_sub)
        
        if mod_index is None:
            print(line)
            raise Exception("Model not found.")
        
        terminals = parts[:mod_index]
        model = parts[mod_index]
        params = parts[(mod_index+1):]

        # Split params and check for params without =
        psplit = []
        for p in params:
            split = p.split("=")
            if len(split)!=2:
                print(line)
                raise Exception("Malformed parameter '"+p+"'.")
            psplit.append(split)
        
        txt = lws + name + " (" + (" ".join(terminals))+") "

        if len(psplit)>0:
            fmted, need_split, split = self.format_params(psplit, len(txt))
            if need_split or split:
                fmted = self.indent(fmted, len(lws)+2)
                txt += "(" + eol + "\n" + fmted
                txt += "\n" + lws + ")"
            else:
                txt += " " + fmted
        
        return txt
