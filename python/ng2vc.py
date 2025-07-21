#!/usr/bin/python3

# Ngspice to VACASK netlist converter
# Use as module, i.e. python3 -m ng2vc ...
# or edit hashbang line and make it executable. 

import sys
from pprint import pprint
from ng2vclib.converter import Converter

if __name__ == "__main__":
    help="""Ngspice to VACASK netlist converter.")
Usage: python3 -m ng2vc [<args>] <input file> [<output file>]
Arguments:
  -h --help  print help
"""
    ndx = 1
    fromFile = None
    toFile = None
    while ndx<len(sys.argv):
        arg = sys.argv[ndx]
        if arg[0]=="-":
            if arg == "--help" or arg == "-h":
                # Print help and exit
                print(help)
                sys.exit(0)
            else:
                print("Unknown argument:", arg)
                print(help)
                sys.exit(1)
        else:
            # Need 1 or 2 more args
            if ndx+1>len(sys.argv):
                print("Need input file.")
                print(help)
                sys.exit(1)
            
            fromFile = arg

            if ndx+2<len(sys.argv):
                print("Too many arguments.")
                print(help)
                sys.exit(1)
            
            if ndx+2==len(sys.argv):
                toFile = sys.argv[ndx+1]
            else:
                toFile = None
            break
    
    cfg = {
        "sourcepath": [ ".", "/home/arpadb/sim/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models" ], 
        "type_map": {
            # name     params             family  remove level  remove version
            "r":     ( {},                "r",     False,        False), 
            "res":   ( {},                "r",     False,        False), 
            "c":     ( {},                "c",     False,        False), 
            "l":     ( {},                "l",     False,        False), 
            "d":     ( {},                "d",     False,        False), 
            "npn":   ( { "type": "1" },   "bjt",   True,         False), 
            "pnp":   ( { "type": "-1" },  "bjt",   True,         False), 
            "njf":   ( { "type": "1" },   "jfet",  True,         False), 
            "pjf":   ( { "type": "-1" },  "jfet",  True,         False), 
            "nmf":   ( { "type": "1" },   "mes",   True,         False), 
            "pmf":   ( { "type": "-1" },  "mes",   True,         False), 
            "nhfet": ( { "type": "1" },   "hemt",  True,         False), 
            "phfet": ( { "type": "-1" },  "hemt",  True,         False), 
            "nmos":  ( { "type": "1" },   "mos",   True,         True), 
            "pmos":  ( { "type": "-1" },  "mos",   True,         True), 
            "nsoi":  ( { "type": "1" },   "soi",   True,         True), 
            "psoi":  ( { "type": "-1" },  "soi",   True,         True), 
        }, 
        "family_map": {
            # family, level, version     file                    module          params
            ("r",     None,  None):    ( "spice/resistor.osdi",  "sp_resistor",  {} ), 
            ("c",     None,  None):    ( "spice/capacitor.osdi", "sp_capacitor", {} ), 
            ("l",     None,  None):    ( "spice/inductor.osdi",  "sp_inductor",  {} ), 
            
            ("d",     None, None):     ( "spice/diode.osdi",     "sp_diode",     {} ), 
            ("d",     1,    None):     ( "spice/diode.osdi",     "sp_diode",     {} ), 
            ("d",     3,    None):     ( "spice/diode.osdi",     "sp_diode",     {} ), 

            ("bjt",   None, None):     ( "spice/bjt.osdi",       "sp_bjt",       {} ), 
            
            ("jfet",  None, None):     ( "spice/jfet1.osdi",     "sp_jfet1",     {} ), 
            ("jfet",  1,    None):     ( "spice/jfet1.osdi",     "sp_jfet1",     {} ), 
            ("jfet",  2,    None):     ( "spice/jfet2.osdi",     "sp_jfet2",     {} ), 
            
            ("mes",   None, None):     ( "spice/mes1.osdi",      "sp_mes1",      {} ), 
            ("mes",   1,    None):     ( "spice/mes1.osdi",      "sp_mes1",      {} ), 

            ("mos",   None, None):     ( "spice/mos1.osdi",      "sp_mos1",      {} ), 
            ("mos",   1,    None):     ( "spice/mos1.osdi",      "sp_mos1",      {} ), 
            ("mos",   2,    None):     ( "spice/mos2.osdi",      "sp_mos2",      {} ), 
            ("mos",   3,    None):     ( "spice/mos3.osdi",      "sp_mos3",      {} ), 
            ("mos",   6,    None):     ( "spice/mos6.osdi",      "sp_mos6",      {} ), 
            ("mos",   9,    None):     ( "spice/mos9.osdi",      "sp_mos9",      {} ), 
            
            ("mos",   8,    None):     ( "spice/bsim3v30.osdi", "sp_bsim3v3",    {} ), 
            ("mos",   8,    "3.3"):    ( "spice/bsim3v30.osdi", "sp_bsim3v3",    {} ), 
            ("mos",   8,    "3.3.0"):  ( "spice/bsim3v30.osdi", "sp_bsim3v3",    {} ), 
            ("mos",   8,    "3.2"):    ( "spice/bsim3v2.osdi",  "sp_bsim3v2",    { "version": '"3.2"'   } ), 
            ("mos",   8,    "3.20"):   ( "spice/bsim3v2.osdi",  "sp_bsim3v2",    { "version": '"3.20"'  } ), 
            ("mos",   8,    "3.2.2"):  ( "spice/bsim3v2.osdi",  "sp_bsim3v2",    { "version": '"3.2.2"' } ), 
            ("mos",   8,    "3.22"):   ( "spice/bsim3v2.osdi",  "sp_bsim3v2",    { "version": '"3.22"'  } ), 
            ("mos",   8,    "3.2.3"):  ( "spice/bsim3v2.osdi",  "sp_bsim3v2",    { "version": '"3.2.3"' } ), 
            ("mos",   8,    "3.23"):   ( "spice/bsim3v2.osdi",  "sp_bsim3v2",    { "version": '"3.23"'  } ), 
            ("mos",   8,    "3.2.4"):  ( "spice/bsim3v2.osdi",  "sp_bsim3v2",    { "version": '"3.2.4"' } ), 
            ("mos",   8,    "3.24"):   ( "spice/bsim3v2.osdi",  "sp_bsim3v2",    { "version": '"3.24"'  } ), 
            ("mos",   8,    "3.1"):    ( "spice/bsim3v1.osdi",  "sp_bsim3v1",    {} ), 
            ("mos",   8,    "3.0"):    ( "spice/bsim3v3.osdi",  "sp_bsim3v0",    {} ), 
            
            ("mos",   49,   None):     ( "spice/bsim3v30.osdi",  "sp_bsim3v3",   {} ), 
            ("mos",   49,   "3.3"):    ( "spice/bsim3v30.osdi",  "sp_bsim3v3",   {} ), 
            ("mos",   49,   "3.3.0"):  ( "spice/bsim3v30.osdi",  "sp_bsim3v3",   {} ), 
            ("mos",   49,   "3.2"):    ( "spice/bsim3v2.osdi",   "sp_bsim3v2",   { "version": '"3.2"'   } ), 
            ("mos",   49,   "3.20"):   ( "spice/bsim3v2.osdi",   "sp_bsim3v2",   { "version": '"3.20"'  } ), 
            ("mos",   49,   "3.2.2"):  ( "spice/bsim3v2.osdi",   "sp_bsim3v2",   { "version": '"3.2.2"' } ), 
            ("mos",   49,   "3.22"):   ( "spice/bsim3v2.osdi",   "sp_bsim3v2",   { "version": '"3.22"'  } ), 
            ("mos",   49,   "3.2.3"):  ( "spice/bsim3v2.osdi",   "sp_bsim3v2",   { "version": '"3.2.3"' } ), 
            ("mos",   49,   "3.23"):   ( "spice/bsim3v2.osdi",   "sp_bsim3v2",   { "version": '"3.23"'  } ), 
            ("mos",   49,   "3.2.4"):  ( "spice/bsim3v2.osdi",   "sp_bsim3v2",   { "version": '"3.2.4"' } ), 
            ("mos",   49,   "3.24"):   ( "spice/bsim3v2.osdi",   "sp_bsim3v2",   { "version": '"3.24"'  } ), 
            ("mos",   49,   "3.1"):    ( "spice/bsim3v10.osdi",  "sp_bsim3v1",   {} ), 
            ("mos",   49,   "3.0"):    ( "spice/bsim3v30.osdi",  "sp_bsim3v0",   {} ), 
        }, 
        "default_models": {
            # letter   osdi file           module
            "r": ( "spice/resistor.osdi",  "sp_resistor" ), 
            "c": ( "spice/capacitor.osdi", "sp_capacitor" ), 
            "l": ( "spice/inductor.osdi",  "sp_inductor" ), 
        }, 
        "default_model_prefix": "default_", 
        "as_toplevel": "auto", 
        "columns": 80, 
        "read_depth": None,    # Fully recursive
        "process_depth": None, # Fully recursive
        "output_depth": None,  
        "all_models": False, 
    }
    converter = Converter(cfg)
    deck = converter.read(fromFile)
    
    # for history, line, control_block in traverse(deck):
    #     lnum, lws, l, eol = line
    #     print(lnum, ":", lws+l, eol)
    # 1/0

    converter.collect_masters()

    # TODO: collect used models and default models
    #       so that we can dump only those that are needed
    #       add annotations to each line
    
    # pprint(converter.data["models"])
    # pprint(converter.data["subckts"])
    # sys.exit(0)
    
    for l in converter.vacask_file():
        print(l)

