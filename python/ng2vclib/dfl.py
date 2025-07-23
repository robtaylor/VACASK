
def default_config(): 
    """
    Returns a default configuration. 
    """
    return {
        "sourcepath": [ "." ], 
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
            ("bjt",   1,    None):     ( "spice/bjt.osdi",       "sp_bjt",       {} ), 
            ("bjt",   4,    None):     ( "spice/vbic.osdi",      "sp_vbic",      {} ), 
            ("bjt",   9,    None):     ( "spice/vbic.osdi",      "sp_vbic",      {} ), 
            
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
        "all_models": False, 
        "default_models": {
            # letter   osdi file           module
            "r": ( "spice/resistor.osdi",  "sp_resistor" ), 
            "c": ( "spice/capacitor.osdi", "sp_capacitor" ), 
            "l": ( "spice/inductor.osdi",  "sp_inductor" ), 
        }, 
        "default_model_prefix": "defmod_", 
        "as_toplevel": "auto", 
        "columns": 80, 
        "read_depth": None,    # Fully recursive
        "process_depth": None, # Fully recursive
        "output_depth": None,  # Fully recursive
    }