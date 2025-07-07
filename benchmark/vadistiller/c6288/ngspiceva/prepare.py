import os

def run(settings):
    vafile = os.path.join(settings["source_dir"], "..", "..", "..", "..", "devices", "spice", "bsim3v30.va")
    return settings["run_openvaf"](vafile, "bsim3v30.osdi")
    
