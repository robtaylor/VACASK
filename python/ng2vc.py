#!/usr/bin/python3

# Ngspice to VACASK netlist converter
# Use as module, i.e. python3 -m ng2vc ...
# or edit hashbang line and make it executable. 

import re
import sys

pat_eolcomment = re.compile(r'(\$|;|//)')
pat_leadspace = re.compile(r'\S')
pat_cidotend = re.compile(r'\.end\s*$', re.IGNORECASE)
pat_spaceseq = re.compile(r' +')
pat_spaceseq = re.compile(r'\s*=\s*')
pat_squotepair = re.compile(r"'(.*?)'")
pat_cbracepair = re.compile(r'\{(.*?)\}')

def convert(fromFile, toFile):
    # Read input file, strip CR/LF and leading whitespace
    try:
        with open(fromFile, 'r') as file:
            lines = [line.lstrip().rstrip('\r\n') for line in file]
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
    for l in lines:
        if len(l)==0:
            # Empty line
            nlines.append(("", l, ""))
            continue
        
        # Separate leading whitespace
        match = pat_leadspace.search(l)
        if match:
            lws = l[:match.start()]
            l = l[match.start():]
        else:
            lws = ""

        # Separate trailing eol comment
        match = pat_eolcomment.search(l)
        if match:
            eolc = l[match.start():]
            l = l[:match.start()]
        else:
            eolc = ""

        if l[0] == "*":
            # Line comment
            comments.append(("", l, ""))
        elif l[0]  == "+":
            # TODO: spaces before +
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
            nlines.append((lws, l, eolc))
    
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
        if l.lower()==".control":
            # Start of control block
            inside_control = True
            continue
        elif l.lower()==".endc":
            # End of control block
            inside_control = False
            continue
        
        # Replace tabs with spaces
        l = l.replace('\t', ' ')
        # Replace sequences of spaces with single space
        l = pat_spaceseq.sub(' ', l)
        # Remove sequences of spaces before and after =
        l = pat_spaceseq.sub('=', l)
        # Remove spaces from strings in single quotes, remove single quotes
        l = pat_squotepair.sub(lambda m: m.group(1).replace(' ', ''), l)
        # Remove spaces from strings in curly braces, remove curly braces
        l = pat_cbracepair.sub(lambda m: m.group(1).replace(' ', ''), l)
        
        # Separate control block, keep case
        if inside_control:
            control_block.append(ll)
            continue
        else:
            # Convert to lowercase and store
            nlines.append((lws, l.lower(), eol))
    
    # Lines now holds pairs of the form (line, eol comment)
    lines = nlines

    # Printout
    print("----")
    for l, ec, eol in lines:
        print(l+ec, eol)
    

    
    

    # Classify lines, collect subcircuit definitions and models
    # A line can be
    # - title
    # - comment
    # - instance
    # - model
    # - subcircuit header
    # - subcircuit end
    # - parameter computation
    # - include
    # - lib
    # - control block

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
    
    if fromFile is None or toFile is None:
        print("Must specify input and output file.")
        print(help)
        sys.exit(1)

    convert(fromFile, toFile)
    sys.exit(0)
