#include <filesystem>
#include "filesystem.h"
#include "filestack.h"
#include "sourceloc.h"
#include "common.h"


namespace NAMESPACE {

std::unordered_map<FileStackIndex,FileStack*> FileStack::idToFs;
FileStackIndex FileStack::at = 1;    // Skip badFileStackId
FileStackIndex FileStack::count = 0;

const FileStackFileIndex FileStack::badFileId = std::numeric_limits<FileStackFileIndex>::max();

const std::string FileStack::blankString = "";

FileStack::FileStack()
    : stringsPool(stringPoolBlockSize, stringPoolGrowthFactor, stringPoolRetries) {
    if (at==maxFileStackId) {
        throw std::length_error("Too many file stacks.");
    }
    id_ = at;
    at++;
    count++;
    idToFs[id_] = this;
}

FileStack::~FileStack() {
    // Remove id to file stack entry
    idToFs.erase(id_);
    // Reduce number of file stacks
    count--;
    // If there are no file stacks in memory, reset id counter
    if (count==0) {
        at = 1;
    }
}

FileStack* FileStack::lookup(FileStackIndex id) {
    auto it = idToFs.find(id);
    if (it!=idToFs.end()) {
        return it->second;
    } else {
        return nullptr;
    }
}

// Look for a file given by fileName
// If the file is specified with absolute path
//   look only at the specified absolute path
// else
//   if the source of the include directive is a file, not a string (obtained from parentId)
//     look in the directory where the parent file is located
//   look in the current working directory
//   look in the include search path
FileStackFileIndex FileStack::addFile(
    const std::string& fileName,  
    const PathList& searchPath, 
    FileStackFileIndex parentId, 
    SourceLineNumber inclusionLine, 
    bool searchParentPath
) {
    // Check limit
    if (stack.size()>=maxFileId) {
        throw std::length_error("Too many files in a file stack.");
    }
    
    // If parent is not given and we have at least one file on the stack, signal an error
    // This means that parsing a file must be done first
    // Later one can only parse strings
    // Maybe remove this limitation in the future
    if (parentId==badFileId && stack.size()>0)
        return badFileId;

    // Absolute paths are tried directly, nothing else is tried
    if (std::filesystem::path(fileName).is_absolute()) {
        std::string canonical;
        if (findFile(fileName, canonical)) {
            stack.emplace_back(
                fileName, 
                "", 
                canonical, 
                nullptr, 
                parentId, 
                inclusionLine
            );
            return stack.size()-1;
        }
    } else {
        // Relative path
        
        std::string canonical;
        
        // First try in parent path
        if (searchParentPath && parentId!=badFileId) {
            // Have parent, get its directory
            std::string parentDirectory = std::move(std::filesystem::path(stack[parentId].canonicalName).parent_path().string());
            if (findFile(fileName, canonical, parentDirectory)) {
                stack.emplace_back(
                    fileName, 
                    "", 
                    canonical, 
                    nullptr, 
                    parentId, 
                    inclusionLine
                );
                return stack.size()-1;
            }
        }
        
        // Try current directory (do not specify path)
        if (findFile(fileName, canonical, "")) {
            stack.emplace_back(
                fileName, 
                "", 
                canonical, 
                nullptr, 
                parentId, 
                inclusionLine
            );
            return stack.size()-1;
        }

        // Finally, try search path
        if (findFile(fileName, canonical, searchPath)) {
            stack.emplace_back(
                fileName, 
                "", 
                canonical, 
                nullptr, 
                parentId, 
                inclusionLine
            );
            return stack.size()-1;
        }
    }
    return badFileId;
}

FileStackFileIndex FileStack::addStringFile(const char* str, size_t n) {
    // Check limit
    if (stack.size()>=maxFileId) {
        throw std::length_error("Too many files in a file stack.");
    }
    
    stack.emplace_back(
        "<string input>", 
        "", 
        "", 
        stringsPool.allocate(str, n), 
        badFileId, 
        0
    );
    return stack.size()-1;
}

FileStackFileIndex FileStack::addStringFile(const char* str) {
    return addStringFile(str, std::strlen(str));
}

FileStackFileIndex FileStack::addStringFile(const std::string& str) {
    return addStringFile(str.c_str(), str.size());
}

FileStackFileIndex FileStack::addStringFile(std::string&& str) {
    return addStringFile(str.c_str(), str.size());
}

FileStackFileIndex FileStack::addRpnString(const char* str, size_t n) {
    // Check limit
    if (rpnStringStack.size()>=maxFileId) {
        throw std::length_error("Too many files in a file stack.");
    }
    
    rpnStringStack.push_back(stringsPool.allocate(str, n));
    return -(rpnStringStack.size()+1);
}

bool FileStack::isString(FileStackFileIndex id) const {
    return id<0;
}

const std::string& FileStack::fileName(FileStackFileIndex id) const {
    if (id>=0) {
        return stack[id].fileName;
    } else {
        return blankString;
    }
}

const std::string& FileStack::sectionName(FileStackFileIndex id) const {
    if (id>=0) {
        return stack[id].sectionName;
    } else {
        return blankString;
    }
}

const std::string& FileStack::canonicalName(FileStackFileIndex id) const {
    if (id>=0) {
        return stack[id].canonicalName;
    } else {
        return blankString;
    }
}

FileStackFileIndex FileStack::parentId(FileStackFileIndex id) const {
    if (id>=0) {
        return stack[id].parentId;
    } else {
        return 0;
    }
}

SourceLineNumber FileStack::inclusionLine(FileStackFileIndex id) const {
    if (id>=0) {
        return stack[id].inclusionLine;
    } else {
        return 0;
    }
}

const char* FileStack::cString(FileStackFileIndex id) const {
    if (id>=0) {
        return stack[id].stringPtr;
    } else {
        return rpnStringStack[-id-1];
    }
}


}
