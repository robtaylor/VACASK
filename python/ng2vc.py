#!/usr/bin/python3

# Ngspice to VACASK netlist converter
# Use as module, i.e. python3 -m ng2vc ...
# or edit hashbang line and make it executable. 

import sys
from pprint import pprint
from ng2vclib.converter import Converter
from ng2vclib import dfl

if __name__ == "__main__":
    help="""Ngspice to VACASK netlist converter.")
Usage: python3 -m ng2vc [<args>] <input file> [<output file>]

If no output file is provided, the converted netlist is printed to 
the standard output. 

Arguments:
  -h --help           print help
  -dr --read-depth    maximal depth to which sources are loaded
                      (infinite by default)
  -dp --process-depth maximal depth to which sources are processed
                      (infinite by default)
  -do --output-depth  maximal depth to which sources are output
                      (infinite ba default)
  -sp --sourcepath    add a directory to the sourcepath where
                      .include and .lib files are found. By default 
                      sourcepath already contains the current directory. 
  -am --all-models    dumps all models, not just those that are used
                      (disabled by default)
"""
    ndx = 1
    fromFile = None
    toFile = None
    read_depth = None
    process_depth = None
    output_depth = None
    sourcepath = [ "." ]
    all_models = False
    while ndx<len(sys.argv):
        arg = sys.argv[ndx]
        if arg[0]=="-":
            if arg == "--help" or arg == "-h":
                # Print help and exit
                print(help)
                sys.exit(0)
            elif arg == "-dr" or arg == "--read-depth":
                if ndx+1>=len(sys.argv):
                    print("Too few arguments.")
                    sys.exit(1)
                ndx += 1
                read_depth = int(sys.argv[ndx])
            elif arg == "-dp" or arg == "--process-depth":
                if ndx+1>=len(sys.argv):
                    print("Too few arguments.")
                    sys.exit(1)
                ndx += 1
                process_depth = int(sys.argv[ndx])
            elif arg == "-do" or arg == "--output-depth":
                if ndx+1>=len(sys.argv):
                    print("Too few arguments.")
                    sys.exit(1)
                ndx += 1
                output_depth = int(sys.argv[ndx])
            elif arg == "-sp" or arg == "--sourcepath":
                if ndx+1>=len(sys.argv):
                    print("Too few arguments.")
                    sys.exit(1)
                ndx += 1
                sourcepath.append(sys.argv[ndx])
            elif arg == "-am" or arg == "--all-models":
                all_models = True
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
        
        ndx += 1
    
    cfg = dfl.default_config()
    cfg["sourcepath"] = [ ".", "/home/arpadb/sim/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models" ]
    cfg["read_depth"] = read_depth
    cfg["process_depth"] = process_depth
    cfg["output_depth"] = output_depth
    cfg["all_models"] = all_models

    converter = Converter(cfg)
    deck = converter.read(fromFile)
    
    # for history, line, control_block in traverse(deck):
    #     lnum, lws, l, eol = line
    #     print(lnum, ":", lws+l, eol)
    # 1/0

    converter.collect_masters()

    # pprint(converter.data["models"])
    # pprint(converter.data["subckts"])
    # sys.exit(0)
    
    out = converter.vacask_file()
    if toFile is None:
        for l in out:
            print(l)
    else:
        with open(toFile, "w") as f:
            for l in out:
                f.write(l)
                f.write("\n")



