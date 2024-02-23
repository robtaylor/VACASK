#ifndef __FILESTACK_DEFINED
#define __FILESTACK_DEFINED

#include <unordered_map>
#include <variant>
#include <vector>
#include <string>
#include <limits>
#include <deque>
#include "pool.h"
#include "common.h"


namespace NAMESPACE {

// File stack holding information on input files thet were parsed 
// in order to construct a circuit. Each FileStack has a unique 
// id. A static repository is used for translating the id to 
// a FileStack object. This is used for reducing the size of Loc 
// objects. 

class FileStack;

typedef struct FileStackEntry {
    std::string fileName;
    std::string sectionName;
    std::string canonicalName;
    const char* stringPtr;
    FileStackFileIndex parentId;
    SourceLineNumber inclusionLine;
    
    FileStackEntry() = delete;
    FileStackEntry(
        const std::string& fileName, 
        const std::string& sectionName, 
        const std::string& canonicalName, 
        const char* stringPtr, 
        FileStackFileIndex parentId, 
        SourceLineNumber inclusionLine) 
        : fileName(fileName), sectionName(sectionName), canonicalName(canonicalName), 
          stringPtr(stringPtr), parentId(parentId), inclusionLine(inclusionLine) {};

    FileStackEntry           (const FileStackEntry&)  = delete;
    FileStackEntry           (      FileStackEntry&&) = default;
    FileStackEntry& operator=(const FileStackEntry&)  = delete;
    FileStackEntry& operator=(      FileStackEntry&&) = default;
} FileStackEntry; 
    
class FileStack {
public:
    typedef std::vector<std::string> PathList;
    static const FileStackIndex maxFileStackId = std::numeric_limits<FileStackIndex>::max();
    static const FileStackIndex badFileStackId = 0;
    static const FileStackFileIndex maxFileId = (std::numeric_limits<FileStackFileIndex>::max()-1)/2;
    static const FileStackFileIndex badFileId;
    
    FileStack();
    ~FileStack();

    FileStack           (const FileStack&)  = delete;
    FileStack           (      FileStack&&) = default;
    FileStack& operator=(const FileStack&)  = delete;
    FileStack& operator=(      FileStack&&) = default;

    // Add a file
    FileStackFileIndex addFile(
        const std::string& fileName,  
        const PathList& searchPath=PathList(), 
        FileStackFileIndex parentId = badFileId, 
        SourceLineNumber inclusionLine=0,
        bool searchParentPath=true
    );

    FileStackIndex id() const { return id_; };

    // Add a string stream, makes sense only for root stack entry
    FileStackFileIndex addStringFile(const char* str, size_t n);
    FileStackFileIndex addStringFile(const char* str);
    FileStackFileIndex addStringFile(const std::string &str);
    FileStackFileIndex addStringFile(std::string &&str);

    // Add a string that is to be parsed as RPN expression
    // Negative indices are reserved for this
    FileStackFileIndex addRpnString(const char* str, size_t n);
    FileStackFileIndex addRpnString(const char* str);
    FileStackFileIndex addRpnString(const std::string &str);
    FileStackFileIndex addRpnString(std::string &&str);

    // Set section name for an entry
    void setSection(FileStackFileIndex id, const std::string& section) { stack[id].sectionName = section; } ;

    // Retrievers
    bool isString(FileStackFileIndex id) const;
    const std::string& fileName(FileStackFileIndex id) const;
    const std::string& sectionName(FileStackFileIndex id) const;
    const std::string& canonicalName(FileStackFileIndex id) const;
    FileStackFileIndex parentId(FileStackFileIndex id) const;
    SourceLineNumber inclusionLine(FileStackFileIndex id) const;
    const char* cString(FileStackFileIndex id) const;

    // Look up FileStack by its id
    static FileStack* lookup(FileStackIndex id);

private:
    // id_ = 0 is reserved for bad FileStack. lookup() resolves it to nullptr. 
    FileStackIndex id_;
    typedef std::vector<FileStackEntry> EntryStack;
    EntryStack stack;

    typedef std::deque<const char*> RpnStringStack;
    RpnStringStack rpnStringStack;

    CStringPool stringsPool;
    
    static std::unordered_map<FileStackIndex,FileStack*> idToFs;
    static FileStackIndex at;
    static FileStackIndex count;

    static const std::string blankString;
};

}

#endif
