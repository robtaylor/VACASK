# Convert IHP PDK for sg13g2 to VACASK format
# For now it does nnot handle the mismatch models, stay tuned. 

# Environmental variables
# PDK_ROOT .. directory created by cloning the PDK 
# PDK .. subdirectory with the PDK, by default ihp-sg13g2

import sys, os, platform, subprocess
from pprint import pprint
from ng2vclib.converter import Converter
from ng2vclib.dfl import default_config

tech_files = [
    # file  read&process depth  output depth
    ( "capacitors_mod.lib", 1, None ), 
    # "capacitors_mod_mismatch.lib", 
    # "capacitors_stat.lib", 
    
    ( "cornerCAP.lib", 0, 0),  
    ( "cornerHBT.lib", 0, 0), 
    ( "cornerMOShv.lib", 0, 0 ), 
    ( "cornerMOSlv.lib", 0, 0 ), 
    ( "cornerRES.lib", 0, 0 ), 
    
    ( "diodes.lib", 1, None ), 
    
    ( "resistors_mod.lib", 1, None ), 
    # "resistors_mod_mismatch.lib", 
    # "resistors_stat.lib", 
    
    ( "sg13g2_bondpad.lib", 1, None ), 
    
    ( "sg13g2_esd.lib", 1, None ), 
    
    ( "sg13g2_hbt_mod.lib", 1, None ), 
    # "sg13g2_hbt_mod_mismatch.lib", 
    # "sg13g2_hbt_stat.lib", 
    
    # "sg13g2_moshv_mismatch.lib", 
    ( "sg13g2_moshv_mod.lib", 1, None ), 
    # "sg13g2_moshv_mod_mismatch.lib", 
    # "sg13g2_moshv_parm.lib", # flattened
    # "sg13g2_moshv_stat.lib", 
    
    # "sg13g2_moslv_mismatch.lib", 
    ( "sg13g2_moslv_mod.lib", 1, None ), 
    # "sg13g2_moslv_mod_mismatch.lib", 
    # "sg13g2_moslv_parm.lib", # flattened
    # "sg13g2_moslv_stat.lib", 
    
    ( "sg13g2_svaricaphv_mod.lib", 1, None),  
    # "sg13g2_svaricaphv_mod_mismatch.lib", 

]

# A bug in Ngspice sg13g2_esd.lib
patches = {
    "sg13g2_esd.lib": [
        (
            ".MODEL diodevss_mod D (tnom = 27 level = 1 is=9.017E-019 rs=200   n=1.03 isr=3.776E-015   ikf=0.0001754 cj0=9.42E-016  m=0.3012  vj=0.6684 bv=11.28 ibv=1E-009 8 nbv=1.324   eg=1.17 xti=3  )", 
            ".MODEL diodevss_mod D (tnom = 27 level = 1 is=9.017E-019 rs=200   n=1.03 isr=3.776E-015   ikf=0.0001754 cj0=9.42E-016  m=0.3012  vj=0.6684 bv=11.28 ibv=1E-009 nbv=1.324   eg=1.17 xti=3  )"
        ), 
    ],
    "sg13g2_svaricaphv_mod.lib": [
        (
            "+ stuac 40", 
            "+ stuac=40"
        )
    ]
}

family_map_update = {
    ("bjt",   4,    None):     ( "vbic_1p3_5t.osdi",     "vbic13_5t",    {} ), 
    ("bjt",   9,    None):     ( "vbic_1p3_5t.osdi",     "vbic13_5t",    {} ), 
}

remove_model_params_update = {
    # translated device name: set of param names
    "vbic13_5t": set(["vbe_max", "vbc_max", "vce_max"])
}

included_va_files = [
    # file                                     options
    ( "psp103/psp103.va",     [ "-D__NGSPICE" ] ), 
    ( "psp103/psp103_nqs.va", [ "-D__NGSPICE" ] ), 
    ( "r3_cmc/r3_cmc.va",     [ "-D__NGSPICE" ] ), 
    ( "mosvar/mosvar.va",     [ "-D__NGSPICE" ] ), 
]

if __name__=="__main__":
    # Get PDK_ROOT environmental variable
    pdkroot = os.getenv("PDK_ROOT")
    if pdkroot is None:
        print("The PDK_ROOT environmental variable must point to the PDK directory.")
        sys.exit(1)

    pdk = os.getenv("PDK")
    if pdk is None:
        pdk = "ihp-sg13g2"

    # Source directory (tech)
    src = os.path.join(pdkroot, pdk, "libs.tech", "ngspice", "models")
    
    # Destination directory (tech)
    dest = os.path.join(pdkroot, pdk, "libs.tech", "vacask", "models")

    # Create destination
    os.makedirs(dest, exist_ok=True)

    # Go through tech files and convert
    osdi_files = set()
    dflmods = set()
    for file, read_process_depth, output_depth in tech_files:
        print("Converting:", file)
        cfg = default_config()
        cfg.update({
            "default_model_prefix": "sg13g2_default_mod_", 
            "sourcepath": [ ".", src ], 
            "read_depth": read_process_depth, 
            "process_depth": read_process_depth, 
            "output_depth": output_depth, 
            "patch": patches, 
        })
        cfg["family_map"].update(family_map_update)
        cfg["remove_model_params"].update(remove_model_params_update)

        srcfile = os.path.join(src, file)
        destfile = os.path.join(dest, file)

        cvt = Converter(cfg)
        cvt.convert(srcfile, destfile)
        
        # OSDI files based on defined and used models
        for mname, in_sub in cvt.data["model_usage"]:
            builtin, mtype, family, level, version, _ = cvt.data["models"][in_sub][mname]
            k = family, level, version
            if k in cvt.cfg["family_map"]:
                file, _, _ = cvt.cfg["family_map"][k]
                osdi_files.add(file)
        
        # OSDI files based on builtin models
        for mt in cvt.data["default_models_needed"]:
            file, module = cvt.cfg["default_models"][mt]
            osdi_files.add(file)
            dflmods.add((mt, module))
    
    # Create .vacaskrc.toml
    print("Creating sample .vacaskrc.toml")
    vacask_cfg="""# VACASK configuration file 
[Paths]
include_path_prefix = [ "$(PDK_ROOT)/$(PDK)/libs.tech/vacask/models" ]
module_path_prefix = [ "$(PDK_ROOT)/$(PDK)/libs.tech/vacask/osdi" ]
"""
    with open(os.path.join(pdkroot, pdk, "libs.tech", "vacask", ".vacaskrc.toml"), "w") as f:
        f.write(vacask_cfg)

    # Platform
    system = platform.system()
    if system=="Windows":
        openvaf_bin = "openvaf-r.exe"
    else:
        openvaf_bin = "openvaf-r"

    # Find this module
    module_path = os.path.dirname(os.path.abspath(__file__))

    # Find OpenVAF
    openvaf = None
    candidates = [
        # <rootdir>/VACASK/python (VS Code build system)
        os.path.join("..", "..", "build.VACASK", "Release", "simulator"), 
        os.path.join("..", "..", "build.VACASK", "Debug", "simulator"), 
        # <rootdir>/lib/vacask/python (Linux)
        os.path.join("..", "..", "..", "bin"), 
        # <rootdir>/lib/python (Windows)
        os.path.join("..", "..", "bin"), 
    ]

    print("Looking for OpenVAF candidate")
    for cand in candidates:
        d = os.path.join(module_path, cand)
        f = os.path.join(d, openvaf_bin)
        print("  ", f)
        if (os.path.isdir(d) and os.path.isfile(f)):
            openvaf = f
            print("Found.")
            break
    
    if openvaf is None:
        print("OpenVAF reloaded not found.")
        sys.exit(1)

    # Modules directory
    mdir = os.path.join(pdkroot, pdk, "libs.tech", "vacask", "osdi")
    os.makedirs(mdir, exist_ok=True)
    
    # Compile modules
    d = os.path.join(pdkroot, pdk, "libs.tech", "verilog-a")
    lead_path = os.path.normpath(mdir)
    for fi, extra_opts in included_va_files:
        f = os.path.join(d, fi)
        fb = os.path.basename(f)
        fo = os.path.join(mdir, fb[:-3]+".osdi")
        
        # Remove leading part from fo, add to list
        fo_n = os.path.normpath(fo)
        fo_r = os.path.relpath(fo_n, lead_path)
        osdi_files.add(fo_r)
        
        print("Compiling", f)
        cmdline = [ openvaf ] + extra_opts + [ "-o", fo, f ]
        retval = subprocess.run(cmdline)
        if retval.returncode != 0:
            print("Verilog-A compiler error.")
            sys.exit(1)

    # Create an include file with common loads and models
    print("Creating common include file")
    txt = ""
    if len(osdi_files)>0:
        txt += "// OSDI files\n"
        for f in osdi_files:
            txt += "load \""+f+"\"\n"
        if len(dflmods)>0:
            txt +="\n"
    if len(dflmods)>0:
        txt += "\n// Default models\n"
        for mt, module in dflmods:
            txt += "model "+cvt.cfg["default_model_prefix"]+mt+" "+module+"\n"
    
    with open(os.path.join(dest, "sg13g2_vacask_common.lib"), "w") as f:
        f.write(txt)



