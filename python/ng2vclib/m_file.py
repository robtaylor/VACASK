import re
import sys
import os
from pprint import pprint

pat_eolcomment = re.compile(r'(\$|;|//)')
pat_leadspace = re.compile(r'^\s+')
pat_cidotend = re.compile(r'\.end\s*$', re.IGNORECASE)
pat_singlequotes = re.compile(r"'(.*?)'")
pat_spaceseq = re.compile(r' +')
pat_spaceequal = re.compile(r'\s*=\s*')
pat_cbracepair = re.compile(r'\{(.*?)\}')

def remove_paired_parentheses_spaces(l):
    """
    Remove spaces from within paired parentheses. 
    """
    stack = []
    result = ""
    for i, char in enumerate(l):
        if char == '(':
            stack.append(i)
            continue
        elif char == ')':
            start = stack.pop()
            if len(stack)==0:
                inner = l[start:i].replace(" ", "")
                result += inner + char
            continue
        
        if len(stack)==0:
            result += char
            continue

    return result

class FileLoaderMixin:
    def preprocess_line(self, line):
        """
        Preprocesses netlist line, ckeck for toplevel circuit. 
        """
        # Convert core line to lowercase
        l = line.lower()

        # Preprocess netlist part (outside control block)
        
        # Strip trailing whitespace from core line
        l = l.rstrip()
        # Replace tabs with spaces
        l = l.replace('\t', ' ')
        # Replace sequences of spaces with single space
        l = pat_spaceseq.sub(' ', l)
        # Remove sequences of spaces before and after =
        l = pat_spaceequal.sub('=', l)
        # Replace pairs of single quotes with curly braces
        l = pat_singlequotes.sub(r'{\1}', l)
        # Remove spaces from strings in curly braces, remove curly braces
        l = pat_cbracepair.sub(lambda m: '{' + m.group(1).replace(' ', '') + '}', l)
        # Remove spaces from within paired parentheses
        l = remove_paired_parentheses_spaces(l)

        # Check if this is a toplevel circuit
        if self.cfg["as_toplevel"]=="auto":
            if l.strip()==".end":
                self.data["is_toplevel"] = True

        return l

    def find_file(self, filename):
        """
        Looks for a file in current directory and in the sourcepath.
        If found, returns its path. 
        If not found, returns None. 
        """
        if os.path.isfile(filename):
            return filename
        
        if not os.path.isabs(filename):
            for p in self.cfg["sourcepath"]:
                full_path = os.path.join(p, filename)
                if os.path.isfile(full_path):
                    return full_path
        
        return None

    def read_file(self, filename, section=None, depth=0, inside_control=False):
        """
        Reads a file. 
        
        If *section* is not None returns only the given section. 
        If section is not found, raises an exception. 
        *depth* is the inclusion deptf of the file. 
        *inside_control* specifies whether the file is inside a control block. 
        
        Returns a tuple with two members if depth>0:
        * *inside_control* status after the end of file is reached
        * deck

        Otherwise returns only the deck. 

        Deck is a tuple with the following members
        * file name
        * section name (for .lib) or None (for .include)
        * a list of tuples (corresponding to lines) with the following members
          * line number
          * leading whitespace
          * core line
          * eol comment

        If recursive_read is enabled in the configuration and a line is an 
        .include or a .lib directive the eol comment entry is a deck. 

        Preprocesses lines, removes trailing whitespace from core line. 
        Merges continued lines (+). 
        Converts all parts to lowercase, except for the control block and 
        .lib/.include directives.  
        Converts pairs of single quotes around expressions into curly braces. 
        Converts tabs to spaces. Removes excessive spaces. 
        Separates control block. 
        """
        # Initialize toplevel circuit status
        self.data["is_toplevel"] = self.cfg["as_toplevel"]=="yes"

        # Read input file, strip CR/LF
        fp = self.find_file(filename)
        if fp is None:
            raise Exception("File "+filename+" not found")
        try:
            with open(fp, 'r') as file:
                lines = [line.rstrip('\r\n') for line in file]
        except:
            raise Exception("Failed to open "+fp)
        
        if depth==0:
            self.data["title"] = lines[0]

        # Split into leading spaces, core, and trailing eol comment
        # Merge lines that start with continuation character. 
        # Line comments between merged lines are dropped. 
        # End of line comment on a line to which a line is merged are dropped. 
        nlines = []
        comments = []
        is_toplevel = False
        collect = section is None
        for lnum, l in enumerate(lines):
            if len(l)==0:
                # Empty line, add to comments
                if collect:
                    comments.append((lnum, "", l, ""))
                continue
            
            # Separate leading whitespace
            match = pat_leadspace.search(l)
            if match:
                lws = l[:match.end()]
                l = l[match.end():]
            else:
                lws = ""

            # Line comment, keep case
            if len(l)>0 and l[0] == "*":
                if collect:
                    comments.append((lnum, lws, l, ""))
                continue
            
            # To lowercase
            llow = l.strip().lower()

            # Is it a .lib or .include line
            islib = llow.startswith(".lib")
            isinclude = llow.startswith(".include")
            isendl = llow.startswith(".endl")
            iscontrol = llow.startswith(".control")
            isendc = llow.startswith(".endc")

            if islib:
                # Check for lib section start
                # Convert sequences of whitespace to single spaces, strip, and split
                lcomp = pat_spaceseq.sub(' ', l.strip()).split(" ")
                # Is it a section marker, i.e. 2 components
                # Section names do not contain whitespace
                if len(lcomp)==2:
                    # Check section
                    sec = lcomp[1]
                    if sec==section:
                        collect = True
                        # Do not store
                        continue
            elif isendl:
                # Check for lib section end
                # Do not check section name, assume file is sane
                if section is not None and collect:
                    collect = False
                    # Do not store
                    continue
            elif iscontrol:
                # Start of control block
                inside_control = True
            elif isendc:
                # End of control block
                inside_control = False
            
            # Not in collect mode, continue
            if not collect:
                continue
            
            # Assume no eol comment
            eolc = ""

            # Recursively read included files
            
            if isinclude:
                if self.cfg.get("recursive_read", False):
                    # Extract file name
                    # Skip .include, strip whitespace
                    s = l[8:].strip()
                    # Unquote
                    if (
                        s.startswith(("'")) and s.endswith(("'")) or 
                        s.startswith(('"')) and s.endswith(('"'))
                    ):
                        s = s[1:-1]
                    # Load include file
                    inside_control, eolc = self.read_file(s, depth=depth+1, inside_control=inside_control)
            elif islib:
                if self.cfg.get("recursive_read", False):
                    # Extract file name and section
                    # Skip .lib, strip whitespace
                    s = l[4:].strip()
                    # Find first whitespace from right to left
                    sndx = s.rfind(" ")
                    if sndx<0:
                        # Not found, must be a stray section marker
                        raise Exception(filename+", line "+str(lnum+1)+": stray section marker.")
                    # Extract file name, strip whitespace
                    lfname = s[:sndx].strip()
                    # Unquote
                    if (
                        lfname.startswith(("'")) and lfname.endswith(("'")) or 
                        lfname.startswith(('"')) and lfname.endswith(('"'))
                    ):
                        lfname = lfname[1:-1]
                    # Get section name
                    sname = s[sndx:].strip()
                    # Load lib file
                    inside_control, eolc = self.read_file(lfname, section=sname, depth=depth+1, inside_control=inside_control)
            elif not inside_control:
                # Not .lib, .include, or control block
                # Separate trailing eol comment
                match = pat_eolcomment.search(l)
                if match:
                    eolc = l[match.start():]
                    l = l[:match.start()]
            
            # Conclude processing
            if len(l)==0:
                # Empty core line, add to comments
                if collect:
                    comments.append((lnum, lws, l, eolc))
                continue
            elif l[0]  == "+":
                # Continuation line
                # Do we have a previous line
                if len(nlines)<=0:
                    raise Exception(filename+", line "+str(lnum+1)+": cannot continue a line without a previous line.")
                # Drop collected comments
                comments = []
                if collect:
                    # Remove end-of-line comment from last line, merge with continuation line and its eol comment
                    prev_lnum, prev_lws, prev_line, _ = nlines[-1]
                    # Add extra space
                    nlines[-1] = (prev_lnum, prev_lws, self.preprocess_line(prev_line+" "+l[1:]), eolc)
                    have_line = True
            else:
                # Ordinary line, flush comments
                if collect:
                    nlines.extend(comments)
                comments = []
                nlines.append((lnum, lws, self.preprocess_line(l), eolc))

        # Flush trailing comments
        nlines.extend(comments)

        # Lines is now a list of tuples of the form 
        # (leading whitespace, core line, trailing eol comment)
        lines = nlines

        if depth>0:
            return (inside_control, (filename, section, lines))
        else:
            if inside_control:
                raise Exception("Unterminated control block found.")

            return (filename, section, lines)
        