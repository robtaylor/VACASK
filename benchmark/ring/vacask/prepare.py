import os

def run(settings):
    vafile = os.path.join(settings["source_dir"], "..", "..", "psp103v4", "vacode", "psp103.va")
    openvaf = settings["run_openvaf"]
    settings["pre_options"] = [ "--skip-embed", "--skip-postprocess", "--no-output" ]
    return openvaf(vafile, "psp103.osdi")
