# Xschem symbol patcher for VACASK

import os, sys, shutil

def simple_patcher(line):
    """
    @pinlist -> ( @pinlist ) 
    """
    if "format" not in line:
        return None
    
    return line.replace("@pinlist", "( @pinlist )").replace("format=", "spectre_format=")

def convert(fname, cvt=simple_patcher):
    """
    If *cvt* is a string, uses exactly that string for spectre_format. 

    If *cvt* is a function, calls it with format line as argument and 
    uses the return value as the spectre_format line. If None is returned 
    conversion fails. 
    """
    # Build orig file name
    origfile = fname+".orig"

    # Does it exists
    if not os.path.isfile(origfile):
        # Not, copy original to .orig
        shutil.copy(fname, origfile)
    
    # Read orig file
    orig_lines = []
    format_line = None
    format_spectre_line = None
    with open(origfile, "r") as f:
        for l in f:
            # Look for format= and format_spectre=
            if l.startswith("format="):
                format_line = len(orig_lines)
            elif l.startswith("spectre_format="):
                format_spectre_line = len(orig_lines)
            
            orig_lines.append(l)

    # No explicitly given format= found, nothing to patch
    if format_line is None:
        os.remove(origfile)
        return

    if isinstance(cvt, str):
        # Explicitly given format
        fstxt = "spectre_format=\""+cvt+"\""
    elif callable(cvt):
        # Convert automatically
        ftxt = orig_lines[format_line]

        # Wrap @pinlist in parentheses
        fstxt = cvt(ftxt)
        if fstxt is None:
            raise Exception("Conversion failed.")
        # fstxt = ftxt.replace("@pinlist", "( @pinlist )").replace("format=", "spectre_format=")
    else:
        # Don't know how to use cvt
        os.remove(origfile)
        raise Exception("Don't know how to use specified converter on "+fname)
        
    # Write format_spectre
    if format_spectre_line is not None:
        # Nothing to do
        pass
    else:
        orig_lines.insert(format_line+1, fstxt)
    
    # Write
    with open(fname, "w") as f:
        for l in orig_lines:
            f.write(l)

if __name__=="__main__":
    pdkroot = os.getenv("PDK_ROOT")
    if pdkroot is None:
        print("The PDK_ROOT environmental variable must point to the PDK directory.")
        sys.exit(1)

    pdk = os.getenv("PDK")
    if pdk is None:
        pdk = "ihp-sg13g2"

    path_pfx = os.path.join(pdkroot, pdk, "libs.tech", "xschem")

    fname = os.path.join(path_pfx, "sg13g2_pr/sg13_hv_nmos.sym")
    convert(fname, None)
    
    