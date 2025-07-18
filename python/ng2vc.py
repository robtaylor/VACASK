#!/usr/bin/python3

# Ngspice to VACASK netlist converter
# Use as module, i.e. python3 -m ng2vc ...
# or edit hashbang line and make it executable. 

import re
import sys
from pprint import pprint
from ng2vclib.converter import Converter
from ng2vclib.generators import traverse, format_history

pat_eolcomment = re.compile(r'(\$|;|//)')
pat_leadspace = re.compile(r'^\s+')
pat_cidotend = re.compile(r'\.end\s*$', re.IGNORECASE)
pat_spaceseq = re.compile(r' +')
pat_spaceequal = re.compile(r'\s*=\s*')
pat_squotepair = re.compile(r"'(.*?)'")
pat_cbracepair = re.compile(r'\{(.*?)\}')

def convert(fromFile, toFile):
    # Read input file, strip CR/LF
    try:
        with open(fromFile, 'r') as file:
            lines = [line.rstrip('\r\n') for line in file]
    except:
        print("Failed to open file", fromFile)
        sys.exit(1)
    
    # Split into leading spaces, core, and trailing eol comment
    # Merge lines that start with continuation character. 
    # Line comments between merged lines are dropped. 
    # End of line comment on a line to which a line is merged are dropped. 
    nlines = []
    comments = []
    is_toplevel = False
    for lnum, l in enumerate(lines):
        if len(l)==0:
            # Empty line, add to comments
            comments.append(("", l, ""))
            continue
        
        # Separate leading whitespace
        match = pat_leadspace.search(l)
        if match:
            lws = l[:match.end()]
            l = l[match.end():]
        else:
            lws = ""

        # Line comment
        if len(l)>0 and l[0] == "*":
            comments.append((lws, l, ""))
            continue
        
        # Separate trailing eol comment
        match = pat_eolcomment.search(l)
        if match:
            eolc = l[match.start():]
            l = l[:match.start()]
        else:
            eolc = ""

        if len(l)==0:
            # Empty line, add to comments
            comments.append((lws, l, eolc))
            continue
        elif l[0]  == "+":
            # Continuation line
            # Do we have a previous line
            if len(nlines)<=0:
                printf("Cannot continue a line without previous line.")
                sys.exit(1)
            # Drop comments
            comments = []
            # Remove end-of-line comment from last line, merge with continuation line and its eol comment
            prev_lws, prev_line, _ = nlines[-1]
            nlines[-1] = (prev_lws, prev_line+l[1:], eolc)
        else:
            # Ordinary line, flush comments
            nlines.extend(comments)
            comments = []
            nlines.append((lws, l, eolc))
    
    # Flush trailing comments
    nlines.extend(comments)

    # Lines is now a list of tuples of the form 
    # (leading whitespace, core line, trailing eol comment)
    lines = nlines

    # Is it a toplevel file (has .end at the end of file)
    for _, l, _ in reversed(lines):
        match = pat_cidotend.search(l)
        if match:
            is_toplevel = True
            break
    
    # Preprocess and extract control block
    nlines = []
    inside_control = False
    control_block = []
    models = { None: [] } # manually specify, if models are define in an included file
    in_sub = None
    
    default_lv = {
        # family  level version
        "r":    ( None, None ), 
        "c":    ( None, None ), 
        "l":    ( None, None ), 
        "d":    ( 1,    None ), 
        "bjt":  ( 1,    None ), 
        "jfet": ( 1,    None ), 
        "mes":  ( 1,    None ), 
        "mos":  ( 1,    None ), 
    }
    
    version_handling = {
        # family  level version      params            version type (None = do not pass)
        # Missing entry .. no special parameters version
        ("d",     None, None):     ( {},               None ), 
        ("d",     1,    None):     ( { "level": "1" }, None ), 
        ("d",     3,    None):     ( { "level": "3" }, None ), 
        
        ("mos",   8,    None):     ( {},               "str" ), 
        ("mos",   8,    "3.3"):    ( {},               "str" ), 
        ("mos",   8,    "3.2"):    ( {},               "str" ), 
        ("mos",   8,    "3.1"):    ( {},               "real" ), 
        ("mos",   8,    "3.0"):    ( {},               None ), 

        ("mos",   49,   None):     ( {},               "str" ), 
        ("mos",   49,   "3.3"):    ( {},               "str" ), 
        ("mos",   49,   "3.2"):    ( {},               "str" ), 
        ("mos",   49,   "3.1"):    ( {},               "real" ), 
        ("mos",   49,   "3.0"):    ( {},               None ), 
    }
    for ll in lines:
        lws, l, eol = ll
        if len(l)==0:
            # Empty line
            nlines.append(ll)
            continue
        elif l[0] == "*":
            # Comment
            nlines.append(ll)
            continue
        
        # Strip trailing whitespace from core line
        l = l.rstrip()
        
        # Check for .control and .endc
        if l.lower().startswith(".control"):
            # Start of control block
            inside_control = True
            continue
        elif l.lower().startswith(".endc"):
            # End of control block
            inside_control = False
            continue
        
        # Replace tabs with spaces
        l = l.replace('\t', ' ')
        # Replace sequences of spaces with single space
        l = pat_spaceseq.sub(' ', l)
        # Remove sequences of spaces before and after =
        l = pat_spaceequal.sub('=', l)
        # Remove spaces from strings in single quotes, remove single quotes
        l = pat_squotepair.sub(lambda m: m.group(1).replace(' ', ''), l)
        # Remove spaces from strings in curly braces, remove curly braces
        l = pat_cbracepair.sub(lambda m: m.group(1).replace(' ', ''), l)
        # Remove spaces from within paired parentheses
        l = remove_paired_parentheses_spaces(l)
        
        # Separate control block, keep case
        if inside_control:
            control_block.append(ll)
            continue
        
        if l.startswith(".subckt"):
            # Detect subckt
            parts = l.split(" ")
            in_sub = parts[1]
        elif l.startswith(".ends"):
            # Detect end of subckt
            in_sub = None
        elif l.startswith(".model"):
            # Collect models
            parts = l.split(" ")
            if in_sub not in models:
                models[in_sub] = []
            models[in_sub].append(parts[1])

            # Detect device type
            module = parts[2]

            # Extract level and version
            # Treat level as integer, version as string


        # Convert to lowercase and store
        nlines.append((lws, l.lower(), eol))
    
    # Lines now holds pairs of the form (line, eol comment)
    lines = nlines

    # Printout
    print("----")
    for lws, l, eol in lines:
        print(lws+l, eol)

    pprint(models)
    

if __name__ == "__main__":
    help="""Ngspice to VACASK netlist converter.")
Usage: python3 -m ng2vc [<args>] <input file> <output file>
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
            # Need exactly 2 more args
            if ndx+2>len(sys.argv):
                print("Too few arguments.")
                print(help)
                sys.exit(1)
            elif ndx+2<len(sys.argv):
                print("Too many arguments.")
                print(help)
                sys.exit(1)
            
            fromFile = arg
            toFile = sys.argv[ndx+1]
            break
    
    cfg = {
        "sourcepath": [ ".", "/home/arpadb/sim/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models" ], 
        "type_map": {
            # name     params             family  remove level  remove version
            "r":     ( {},                "r",     False,        False), 
            "res":   ( {},                "r",     False,        False), 
            "c":     ( {},                "c",     False,        False), 
            "l":     ( {},                "l",     False,        False), 
            "d":     ( {},                "d",     True,         False), 
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
            # family, level, version     file                    module
            ("r",     None,  None):    ( "spice/resistor.osdi",  "sp_resistor" ), 
            ("c",     None,  None):    ( "spice/capacitor.osdi", "sp_capacitor" ), 
            ("l",     None,  None):    ( "spice/inductor.osdi",  "sp_inductor" ), 
            
            ("d",     None, None):     ( "spice/diode.osdi",     "sp_diode" ), 
            ("d",     1,    None):     ( "spice/diode.osdi",     "sp_diode" ), 
            ("d",     3,    None):     ( "spice/diode.osdi",     "sp_diode" ), 

            ("bjt",   None, None):     ( "spice/bjt.osdi",       "sp_bjt" ), 
            
            ("jfet",  None, None):     ( "spice/jfet1.osdi",     "sp_jfet1" ), 
            ("jfet",  1,    None):     ( "spice/jfet1.osdi",     "sp_jfet1" ), 
            ("jfet",  2,    None):     ( "spice/jfet2.osdi",     "sp_jfet2" ), 
            
            ("mes",   None, None):     ( "spice/mes1.osdi",      "sp_mes1" ), 
            ("mes",   1,    None):     ( "spice/mes1.osdi",      "sp_mes1" ), 

            ("mos",   None, None):     ( "spice/mos1.osdi",      "sp_mos1" ), 
            ("mos",   1,    None):     ( "spice/mos1.osdi",      "sp_mos1" ), 
            ("mos",   2,    None):     ( "spice/mos2.osdi",      "sp_mos2" ), 
            ("mos",   3,    None):     ( "spice/mos3.osdi",      "sp_mos3" ), 
            ("mos",   6,    None):     ( "spice/mos6.osdi",      "sp_mos6" ), 
            ("mos",   9,    None):     ( "spice/mos9.osdi",      "sp_mos9" ), 
            
            ("mos",   8,    None):     ( "spice/bsim3v30.osdi", "sp_bsim3v3" ), 
            ("mos",   8,    "3.3"):    ( "spice/bsim3v30.osdi", "sp_bsim3v3" ), 
            ("mos",   8,    "3.3.0"):  ( "spice/bsim3v30.osdi", "sp_bsim3v3" ), 
            ("mos",   8,    "3.2"):    ( "spice/bsim3v2.osdi",  "sp_bsim3v2" ), 
            ("mos",   8,    "3.2.2"):  ( "spice/bsim3v2.osdi",  "sp_bsim3v2" ), 
            ("mos",   8,    "3.2.3"):  ( "spice/bsim3v2.osdi",  "sp_bsim3v2" ), 
            ("mos",   8,    "3.2.4"):  ( "spice/bsim3v2.osdi",  "sp_bsim3v2" ), 
            ("mos",   8,    "3.1"):    ( "spice/bsim3v1.osdi",  "sp_bsim3v1" ), 
            ("mos",   8,    "3.0"):    ( "spice/bsim3v3.osdi",  "sp_bsim3v0" ), 
            
            ("mos",   49,   None):     ( "spice/bsim3v30.osdi",  "sp_bsim3v3" ), 
            ("mos",   49,   "3.3"):    ( "spice/bsim3v30.osdi",  "sp_bsim3v3" ), 
            ("mos",   49,   "3.3.0"):  ( "spice/bsim3v30.osdi",  "sp_bsim3v3" ), 
            ("mos",   49,   "3.2"):    ( "spice/bsim3v24.osdi",  "sp_bsim3v2" ), 
            ("mos",   49,   "3.2.2"):  ( "spice/bsim3v24.osdi",  "sp_bsim3v2" ), 
            ("mos",   49,   "3.2.3"):  ( "spice/bsim3v24.osdi",  "sp_bsim3v2" ), 
            ("mos",   49,   "3.2.4"):  ( "spice/bsim3v24.osdi",  "sp_bsim3v2" ), 
            ("mos",   49,   "3.1"):    ( "spice/bsim3v10.osdi",  "sp_bsim3v1" ), 
            ("mos",   49,   "3.0"):    ( "spice/bsim3v30.osdi",  "sp_bsim3v0" ), 
        }, 
        "as_toplevel": "auto", 
        "columns": 80, 
        "read_depth": None,    # Fully recursive
        "process_depth": None, # Fully recursive
        "output_depth": 0,     # Toplevel only
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

    sys.exit(0)


    if fromFile is None or toFile is None:
        print("Must specify input and output file.")
        print(help)
        sys.exit(1)

    convert(fromFile, toFile)
    sys.exit(0)
