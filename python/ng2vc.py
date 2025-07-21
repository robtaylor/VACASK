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
    
    converter = Converter()
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
    
    out = converter.vacask_file()
    if toFile is None:
        for l in out:
            print(l)
    else:
        with open(toFile, "w") as f:
            for l in out:
                f.write(l)
                f.write("\n")



