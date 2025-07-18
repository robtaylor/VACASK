
# A generator that traverses a deck
def traverse(deck, recursive=True, input_history=[], parent_line=None, inside_control=False):
    """
    A generator that traverses a *deck*. 
    if *recursive* is true it also traverses included files. 

    *input_history* is the history up to the *deck*. 

    *parent_line* is the line number on which this *deck* was included. 

    *inside_control* is a flag specifying whether the deck is inside 
    a .control block. 

    It yields 
    * a list of history items which are tuples holding
      * line number in parent file where this file was included, 
        None for the toplevel file. 
      * filename
      * section
    * a line which is a tuple holding
      * line number
      * leading whitespace
      * core line
      * eol comment .. a string or a tuple if the line corresponds 
        to a .include or .lib directive
    * a flag indicating this line is a part of control clock
    """
    filename, section, lines = deck
    history = input_history + [ (parent_line, filename, section) ]
    
    for line in lines:
        lnum, lws, l, eol = line
        # Check for control
        ll = l.lower()
        if ll.startswith(".control"):
            inside_control = True
        elif ll.startswith(".endc"):
            inside_control = False

        # Is it a normal line
        if not isinstance(eol, tuple):
            # Yes
            yield (history, line, inside_control)
        else:
            # No, recurse if requested
            if recursive:
                for subh, subl, inside_control in traverse(eol, recursive, history, lnum, inside_control=inside_control):
                    yield (subh, subl, inside_control)

def format_history(history, lineno):
    """
    Formats *history* for an error on line number *lineno*. 
    """
    txt = "On line "+str(lineno+1)+" of "
    for h in reversed(history):
        pline, fname, sec = h
        txt += "'"+fname+"'"
        if pline is not None:
            if sec is not None:
                txt += "included as section"+sec+" on line "+str(pline+1)+" of "
            else:
                txt += "included on line "+str(pline+1)+" of "
    
    return txt
    