#include <exception>
#include <fstream>
#include <ios>
#include "sourceloc.h"
#include "filesystem.h"
#include "common.h"


namespace NAMESPACE {

const Loc Loc::bad = Loc();

Loc::Loc() 
    : entry(Entry(FileStack::badFileStackId, FileStack::badFileId, 0, 0)) {
}

void Loc::addEntry(FileStackIndex fs, FileStackFileIndex file, SourceLineNumber l, SourceColumnNumber o) {
    entry = Entry(fs, file, l, o);
}

Loc::Loc(FileStackIndex fs, FileStackFileIndex file, SourceLineNumber l, SourceColumnNumber o) {
    addEntry(fs, file, l, o);
};

Loc::Loc(const FileStack* fs, FileStackFileIndex file, SourceLineNumber l, SourceColumnNumber o) {
    addEntry(fs->id(), file, l, o);
}

// Loc::Loc(const Position& pos) {
//     addEntry(pos.fileStack->id(), pos.fileStackPosition, pos.line, pos.offset);
// }
// 
// Loc::Loc(const Location& loc) {
//     addEntry(loc.begin.fileStack->id(), loc.begin.fileStackPosition, loc.begin.line, loc.begin.offset);
// }

std::tuple<const FileStack*, FileStackFileIndex, SourceLineNumber, SourceColumnNumber> Loc::data() const { 
    auto fs = FileStack::lookup(entry.fileId_);
    return std::make_tuple(FileStack::lookup(entry.fileStackId_), entry.fileId_, entry.line_, entry.offset_); 
};

std::string Loc::toString() const {
    std::ostringstream s;
    std::string str;

    SourceColumnNumber column;
    auto [fileStack, fileId, line, offset] = data();

    const SourceColumnNumber maxFront = 35;
    const SourceColumnNumber maxLen = 70;

    // Get line
    bool ok;
    if (fileStack->canonicalName(fileId).size()>0) {
        std::ifstream f(fileStack->canonicalName(fileId).c_str());
        auto [found, col] = getLineByOffset(f, offset, str);
        column = col;
        ok = found;
    } else {
        auto [found, col] = getLineByOffset(fileStack->cString(fileId), offset, str);
        column = col;
        ok = found;
    }
    
    if (ok) {
        auto col = column;
        SourceColumnNumber nErase = 0, nEraseBack = 0;
        s << "  ";
        if (col>maxFront) {
            // Trim front 
            nErase = col-maxFront;
            str.erase(0, nErase);
            s << "... ";
            col -= nErase; // length of deleted part
            col += 4; // length of elipsis
        }
        if (str.size()>maxLen) {
            // Trim back
            nEraseBack = str.size()-maxLen;
            str.erase(str.size()-nEraseBack, nEraseBack);
        }
        // Output (trimmed)
        s << str;
        if (nEraseBack) {
            s << " ...";
        }
        s << "\n" << "  ";
        for(SourceColumnNumber i=1; i<col; i++) 
            s << " ";
        s << "^";
    }

    if (entry.fileStackId_!=FileStack::badFileStackId) {
        // <line>:<column> in <filename>[, section <section>]
        // included on line <lineno> of <filename>[, section <section>] 
        // included on line ...
        s << "\n";

        // Is it just a string with no inclusion history
        if (fileStack->isString(fileId)) {
            return std::string("<string input>");
        }
        s << line << ":" << column << " (0x" << std::hex << offset << std::dec << ")" << " in ";
        auto pos = fileId;
        while (pos!=FileStack::badFileId) {
            FileStackFileIndex next = fileStack->parentId(pos);
            
            if (fileStack->canonicalName(pos).size()>0) {
                s << fileStack->canonicalName(pos);
            } else {
                s << fileStack->fileName(pos);
            }
            if (fileStack->sectionName(pos).size()>0) {
                s << ", section " << fileStack->sectionName(pos);
            }
            if (next!=FileStack::badFileId) {
                s << "\n";    
            }
            if (fileStack->inclusionLine(pos)>0) {
                s << "  included on line " << fileStack->inclusionLine(pos) << " of ";
            }
            pos = next;
        }
    }

    return s.str();
}

std::ostream& operator<<(std::ostream& ostr, const Loc& p) {
    ostr << p.toString();
    return ostr;
}

}
