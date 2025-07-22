
class ModelMixin:
    def process_model(self, lws, line, eol, in_sub):
        """
        Process a model line. 
        """
        # Model
        name = line.split(" ")[1]
        builtin, mtype, family, level, version, params = self.data["models"][in_sub][name]

        paren = False
        if mtype in self.cfg["type_map"]:
            # Builtin model
            extra_params, family, _, _ = self.cfg["type_map"][mtype]
            osdifile, module, extra_family_params = self.cfg["family_map"][family, level, version]
            vcline = "model "+name+" "+module
            if len(extra_params)>0 or len(params)>0:
                vcline += " ( "+eol
                paren = True
            else:
                vcline += " "+eol

            if len(extra_params)>0:
                vcline += (
                    "\n"+self.format_extra_params(extra_params, lws, 0)+
                    self.format_extra_params(extra_family_params, " ", 0)
                )
        else:
            vcline = "model "+name+" "+mtype
            if len(params)>0:
                vcline += " ( "+eol
                paren = True
            else:
                vcline += " "+eol
        
        if len(params)>=0:
            fmted, _, _ = self.format_params(params, len(vcline))
            fmted = self.indent(fmted, len(lws)+2)
            vcline += "\n" + fmted
        
        if paren:
            vcline += "\n" + lws + ")"

        return vcline
        
