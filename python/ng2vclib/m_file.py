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
    def find_file(self, filename):
        """
        Looks for a file in current directory and in the sourcepath.
        If found, returns the canonical path. 
        If not found, returns None. 
        """
        if os.path.isfile(filename):
            return os.path.realpath(filename)
        
        if not os.path.isabs(filename):
            for p in self.cfg["sourcepath"]:
                full_path = os.path.join(p, filename)
                if os.path.isfile(full_path):
                    return os.path.realpath(full_path)
        
        return None

    def read_file(self, filename, recursive=True, section=None, depth=0):
        """
        Reads a file. If *recursive* is True, recursively handles 
        .include and .lib directives. 
        
        If *section* is not None returns only the given section. 
        If section is not found, raises an exception. 
        
        Returns a list of tuples corresponding to lines
        with the following members:
        * leading whitespice
        * core line
        * eol comment

        If a line is an .iclude or a .lib directive and *recursive* is 
        True the eol comment entry is a list holding 
        * the included file's lines
        * file name
        * optional section

        Preprocesses lines, removes trailing whitespace from core line. 
        Merges continued lines (+). 
        Converts all parts to lowercase, except for the control block and 
        .lib/.include directives.  
        Converts pairs of single quotes around expressions into curly braces. 
        Converts tabs to spaces. Removes excessive spaces. 
        Separates control block. 
        """
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
                if collect:
                    comments.append((lws, l, ""))
                continue
            
            # Separate trailing eol comment
            match = pat_eolcomment.search(l)
            if match:
                eolc = l[match.start():]
                l = l[:match.start()]
            else:
                eolc = ""

            # Check for lib section start
            llow = l.strip().lower()
            if llow.startswith(".lib"):
                lcomp = pat_spaceseq.sub(' ', l.strip()).split(" ")
                # Is it a section marker
                if len(lcomp)==2:
                    # Check section
                    sec = lcomp[1]
                    if sec==section:
                        collect = True
                        # Do not store
                        continue
            
            # Check for lib section end
            if llow.startswith(".endl"):
                # Do not check section name, assume file is sane
                if section is not None and collect:
                    collect = False
                    # Do not store
                    continue

            if len(l)==0:
                # Empty line, add to comments
                if collect:
                    comments.append((lws, l, eolc))
                continue
            elif l[0]  == "+":
                # Continuation line
                # Do we have a previous line
                if len(nlines)<=0:
                    raise Exception(filename+", line "+str(lnum+1)+": cannot continue a line without previous line.")
                # Drop comments
                comments = []
                if collect:
                    # Remove end-of-line comment from last line, merge with continuation line and its eol comment
                    prev_lws, prev_line, _ = nlines[-1]
                    # Add extra space
                    nlines[-1] = (prev_lws, prev_line+" "+l[1:], eolc)
            else:
                # Ordinary line, flush comments
                if collect:
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
            # Replace pairs of single quotes with curly braces
            l = pat_singlequotes.sub(r'{\1}', l)
            # Remove spaces from strings in curly braces, remove curly braces
            l = pat_cbracepair.sub(lambda m: '{' + m.group(1).replace(' ', '') + '}', l)
            # Remove spaces from within paired parentheses
            l = remove_paired_parentheses_spaces(l)
            
            # Separate control block, keep case
            if inside_control:
                self.data["control"].append(ll)
                continue

            # Convert to lowercase if not .lib or .include
            llow = l.lower()
            if ( 
                not llow.startswith(".include") and
                not llow.startswith(".lib")
            ):
                l = llow
            
            # Recursively read
            if recursive:
                if llow.startswith(".include"):
                    # Extract file name
                    s = l.split(" ")[1]
                    # Unquote - do not use spaces in file names
                    if s.startswith(("'", '"')) and s.endswith(("'", '"')):
                        s = s[1:-1]
                    # Load include file
                    eol = [ self.read_file(s, recursive, depth=depth+1), s ]
                elif llow.startswith(".lib"):
                    # Extract file name and section
                    lcomp = l.split(" ")
                    s1 = lcomp[1]
                    s2 = lcomp[2]
                    # Unquote - do not use spaces in file names
                    if s1.startswith(("'", '"')) and s1.endswith(("'", '"')):
                        s1 = s1[1:-1]
                    if s2.startswith(("'", '"')) and s2.endswith(("'", '"')):
                        s2 = s2[1:-1]
                    # Load lib file
                    eol = [ self.read_file(s1, recursive, section=s2, depth=depth+1), s1, s2 ]

            # Store
            nlines.append((lws, l, eol))
        
        # Lines now holds pairs of the form (line, eol comment)
        lines = nlines

        # Detect toplevel circuit
        self.data["is_toplevel"] = False
        if self.cfg["as_toplevel"]=="auto":
            for _, l, _ in lines:
                if l.strip() == ".end":
                    self.data["is_toplevel"] = True
                    break

        else:
            self.data["is_toplevel"] = self.cfg["as_toplevel"]=="yes"
        

        self.data["lines"] = lines

        return lines        