import os

def run(settings):
    vafile = os.path.join(settings["source_dir"], "..", "..", "psp103v4", "vacode", "psp103.va")
    settings["pre_options"] = [ "--skip-embed", "--skip-postprocess", "--no-output" ]
    return True
