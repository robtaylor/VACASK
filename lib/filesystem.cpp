
#include <fstream>
#include <filesystem>
#include <string>
#include "filesystem.h"
#include "libplatform.h"
#include "common.h"


namespace NAMESPACE {

// Try to find file only if absolute path to file is given
bool findFile(const std::string& fileName, std::string& canonicalPath) {
    try {
        // Absolute path to file given, ignore directory
        if (std::filesystem::path(fileName).is_absolute()) {
            auto canonical = std::filesystem::canonical(fileName);
            if (std::filesystem::exists(canonical)) {
                if (std::filesystem::file_size(canonical)>std::numeric_limits<SourceColumnNumber>::max()) {
                    throw std::length_error(std::string("Input file '")+canonical.string()+"' too big.");
                }
                canonicalPath = std::move(canonical.string());
                return true;
            } 
        } 
        return false;
    } catch (std::filesystem::filesystem_error &fe) {
        return false;
    }
}

bool findFile(const std::string& fileName, std::string& canonicalPath, const std::string& directory) {
    try {
        if (std::filesystem::path(fileName).is_absolute()) {
            // Absolute path to file given, ignore directory
            return findFile(fileName, canonicalPath);
        } else {
            // Relative path, use given directory
            std::filesystem::path childPath;
            if (directory.size()>0) {
                // Path has nonzero length, use it as starting point
                childPath = std::filesystem::canonical(directory);
            } else {
                // Zero length path, use current directory as starting point
                auto childPath = std::filesystem::current_path();
            }
            // Add file name, see if it exists
            childPath /= fileName;
            auto canonical = std::filesystem::canonical(childPath);
            if (std::filesystem::exists(canonical)) {
                if (std::filesystem::file_size(canonical)>std::numeric_limits<SourceColumnNumber>::max()) {
                    throw std::length_error(std::string("Input file '")+canonical.string()+"' too big.");
                }
                canonicalPath = std::move(canonical.string());
                return true;
            }
        }
        return false;
    } catch (std::filesystem::filesystem_error &fe) {
        return false;
    }
}

std::tuple<bool, std::string> findFileInSystemPath(const std::string& fileName) {
    for(auto& dir : LibPlatform::systemPath()) {
        std::string canonicalPath;
        if (findFile(fileName, canonicalPath, dir)) {
            return std::make_tuple(true, canonicalPath);
        }
    }
    return std::make_tuple(false, "");
}

// Try list of directories, also handle absolute path
bool findFile(const std::string& fileName, std::string& canonicalPath, const std::vector<std::string>& path) {
    for(auto it=path.begin(); it!=path.end(); ++it) {
        if (findFile(fileName, canonicalPath, *it))
            return true;
    }
    return false;
}

// No longer used
bool getLine(std::istream& f, int line, std::string &destination) {
    std::string s;

    for (int i=0; i<line-1; i++) {
        f.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        if (f.eof()) 
            return false;
    }
    
    std::getline(f, destination);
    return true;
}

// No longer used
bool getLine(const char* s, int line, std::string &destination) {
    // Find line
    int atLine = 1;
    const char *ptr;
    for(ptr=s; *ptr; ptr++) { 
        if (atLine == line) {
            break;
        }
        if (*ptr=='\n') {
            atLine++;
        }
    }
    
    // Not enough lines
    if (*ptr==0) {
        return false;
    }
    
    // Find end of line
    const char* endPtr;
    for(endPtr=ptr; *endPtr; endPtr++) {
        if (*endPtr=='\n') {
            break;
        }
    }

    destination = std::move(std::string(ptr, endPtr-ptr));
    return true;
}

std::tuple<bool, std::streamoff> getLineByOffset(std::istream& f, std::streamoff offset, std::string& destination) {
    std::streamoff back = 1;
    std::streamoff from = offset;
    std::streamoff to;
    std::streamoff beginning = 0; // Start of file is start of file by default
    bool found = false;
    
    // Find start of line
    while (!found) {
        to = from;
        from = to - back;
        if (from<0) {
            from = 0;
        }

        f.seekg(from);
        if (f.fail()) {
            return std::make_tuple(false, 0);
        }
        for(auto i=from; i<to; i++) {
            auto ch = f.get();
            if (ch=='\n') {
                // Found newline
                found = true;
                beginning = f.tellg();
            }
        }
        if (!found) {
            back *= 2;
        }

        if (to==0) {
            // Reached start of file, give up
            break;
        }
    }

    // Read line
    f.seekg(beginning);
    if (f.fail()) {
        return std::make_tuple(false, 0);
    }
    std::getline(f, destination);

    // First column has number 1
    return std::make_tuple(true, offset-beginning+1);
}

std::tuple<bool, std::streamoff> getLineByOffset(const char* s, std::streamoff offset, std::string& destination) {
    std::streamoff beginning = 0;
    std::streamoff end = 0;
    if (offset>0) {
        // Find start of line
        for(auto i=offset-1; i>=0; i--) {
            if (s[i]=='\n') {
                beginning = i+1;
                break;
            }
            if (i==0) {
                // Counter might be unsigned, stop after 0 is processed
                break;
            }
        }
    }
    // Find end of string
    for(end=offset; s[end]; end++) {
        if (s[end]=='\n') {
            break;
        }
    }
    // Build string
    destination = std::move(std::string(s+beginning, end-beginning));
    return std::make_tuple(true, offset-beginning+1);
}

}
