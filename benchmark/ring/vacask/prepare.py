import os

def run(settings):
    vafile = os.path.join(settings["source_dir"], "..", "..", "psp103v4", "vacode", "psp103.va")
    openvaf = settings["run_openvaf"]
    return openvaf(vafile, "psp103.osdi")
