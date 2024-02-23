#ifndef __FILESYSTEM_DEFINED
#define __FILESYSTEM_DEFINED

#include <iostream>
#include <string>
#include <variant>
#include <vector>
#include "common.h"


namespace NAMESPACE {

// Find file, try only if absolute path is given
bool findFile(const std::string& fileName, std::string& canonicalPath);

// Find file, return its canonical path
// If absolute path is given try that (use previous function), else try to find file in given directory
// Directory name of length 0 corresponds to current directory
bool findFile(const std::string& fileName, std::string& canonicalPath, const std::string& directory);

// Look for a file in a set of paths, use previous function to search one directory
bool findFile(const std::string& fileName, std::string& canonicalPath, const std::vector<std::string>& path);

// Look for a file in system path
std::tuple<bool, std::string> findFileInSystemPath(const std::string& fileName);

// Read a specific line from a stream
bool getLine(std::istream&, int line, std::string &destination);

// Read a specific line from a c string
bool getLine(const char* s, int line, std::string &destination);

// Read a line containing a given offset from a file
// CHecks if offset is valid
// Compute column position, does not work with UTF-8
std::tuple<bool, std::streamoff> getLineByOffset(std::istream&, std::streamoff offset, std::string& destination);

// Read a line containing a given offset from a c string
// Does not check if offset is valid
// Compute column position, does not work with UTF-8
// 
std::tuple<bool, std::streamoff> getLineByOffset(const char* s, std::streamoff offset, std::string& destination);

}

#endif
