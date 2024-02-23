#ifndef __SOURCELOC_DEFINED
#define __SOURCELOC_DEFINED

#include <deque>
#include <unordered_map>
#include <limits>
#include "filestack.h"
#include "common.h"


namespace NAMESPACE {

// Contrary to identifiers, locations are typically referenced only once. 
// They are also unique across multiple circuits. 
// Therefore it would not make sense to have a repository common to all circuits. 
// But there is one reason: Rpn::Entry is much smaller if each Loc is actually 
// a location pointer (8 bytes). We need to keep a separate repository 
// for each FileStack we encounter so we can easily delete Entries
// We replace a FileStack pointer with a FileStackId (uint32_t) which makes 
// the Loc structure smaller (FileStackId+FileStackPosition+Line+Column =
// 4+4+4+4 = 16 bytes in size). If we use a pointer to FileStack, Loc size is 20 
// which (if counting for alignment) effectively results in 24 bytes. 
class Loc {
public:
    typedef struct Entry {
        FileStackIndex fileStackId_;
        FileStackFileIndex fileId_;
        SourceLineNumber line_;
        SourceColumnNumber offset_;
        Entry() : fileStackId_(FileStack::badFileStackId), fileId_(FileStack::badFileId), line_(0), offset_(0) {};
        Entry(FileStackIndex fileStackId, FileStackFileIndex fileId, SourceLineNumber line, SourceColumnNumber offset) 
            : fileStackId_(fileStackId), fileId_(fileId), line_(line), offset_(offset) {};
    } Entry;

    Loc           (const Loc&)  = default;
    Loc           (      Loc&&) = default;
    Loc& operator=(const Loc&)  = default;
    Loc& operator=(      Loc&&) = default;

    // Default constructor creates a bad location object
    Loc();

    // Constructs a Loc object from FileStack id, file id, line, and column
    Loc(FileStackIndex fs, FileStackFileIndex file, SourceLineNumber l, SourceColumnNumber offset);

    // Constructs a Loc object from FileStack ptr, file id, line, and column
    Loc(const FileStack* fs, FileStackFileIndex file, SourceLineNumber l, SourceColumnNumber offset);

    // Constructs a Loc object from a Bison Position object
    // Loc(const Position& pos);

    // Constructs a Loc object from a Bison Location object
    // Loc(const Location& loc);
    
    // Contents of location object
    std::tuple<const FileStack*, FileStackFileIndex, SourceLineNumber, SourceColumnNumber> data() const;

    // Format as string
    std::string toString() const;
    // std::string trace() const;

    friend std::ostream& operator<<(std::ostream& ostr, const Loc& p);

    inline bool operator==(const Loc& other) const {
        return entry.fileStackId_==other.entry.fileStackId_ &&
            entry.fileId_==other.entry.fileId_ &&
            entry.line_==other.entry.line_ &&
            entry.offset_==other.entry.offset_; 
    };
    inline bool operator!=(const Loc& other) const { return !(*this==other); };

    inline operator bool() const { return entry.fileStackId_ != FileStack::badFileStackId; };
    
    static const Loc bad;

private:    
    void addEntry(FileStackIndex fs, FileStackFileIndex file, SourceLineNumber l, SourceColumnNumber offset);
    Entry entry;
};

}

#endif