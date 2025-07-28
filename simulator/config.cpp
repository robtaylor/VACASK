#include "status.h"
#include "simulator.h"
#include "platform.h"

#include <toml++/toml.h>
#include <iostream>
#include <string>
#include <cstdlib>

#include "common.h"


namespace NAMESPACE {

void replaceEnvVariables(std::string& input) {
    size_t pos = 0;
    while ((pos = input.find("$(", pos)) != std::string::npos) {
        size_t endPos = input.find(')', pos);
        if (endPos == std::string::npos) {
            break; // No closing parenthesis
        }

        std::string varName = input.substr(pos + 2, endPos - pos - 2); // Extract NAME from $(NAME)
        const char* envValue = std::getenv(varName.c_str());
        std::string replacement = envValue ? envValue : "";

        input.replace(pos, endPos - pos + 1, replacement);
        // Don't advance pos; the replacement might contain more variables
    }
}

// OK, found, string
std::tuple<bool, bool, std::string> getString(toml::v3::table* db, std::string name, Status& s) {
    if (auto node = (*db)[name]) {
        if (node.is_string()) {
            std::string str = node.as_string()->get();
            replaceEnvVariables(str);
            return std::make_tuple(true, true, str);
            
        } else {
            s.set(Status::Syntax, std::string("Expected string for ")+name+".");
            return std::make_tuple(false, true, "");
        }
    }
    return std::make_tuple(true, false, "");
}

// OK, array ptr
std::tuple<bool, toml::v3::array*> getStringArray(toml::v3::table* db, std::string name, Status& s) {
    auto path = (*db)[name];
    if (path) {
        if (auto arr = path.as_array()) {
            for (size_t i = 0; i < arr->size(); ++i) {
                if (!(*arr)[i].is_string()) {
                    s.set(Status::Syntax, std::string("Non-string entry at index ") + std::to_string(i) + " in " + name+".");
                    return std::make_tuple(false, nullptr);
                }
            }
            return std::make_tuple(true, arr);
        } else {
            // Not an array
            s.set(Status::Syntax, std::string("Key ") + name + " is not an array.");
            return std::make_tuple(false, nullptr);
        }
    }

    return std::make_tuple(true, nullptr);
}

std::vector<std::string> arrayToStringVector(toml::v3::array* arr) {
    std::vector<std::string> vec;
    for (const auto& node : *arr) {
        if (auto val = node.value<std::string>()) {
            replaceEnvVariables(*val);
            vec.push_back(*val);
        }
    }
    return vec;
}

// Return value: 
bool readConfig(std::ifstream& input, const std::string& filename, Status& s) {
    try {
        // Parse TOML
        auto config = toml::parse(input);

        if (auto database = config["Paths"].as_table()) {
            if (auto [ok, arr] = getStringArray(database, "include_path_prefix", s); ok) {
                // Prefix path
                if (arr) {
                    Simulator::prependIncludePath(std::move(arrayToStringVector(arr)));
                }
            } else {
                return false;
            }
            if (auto [ok, arr] = getStringArray(database, "include_path_suffix", s); ok) {
                // Suffix path
                if (arr) {
                    Simulator::appendIncludePath(std::move(arrayToStringVector(arr)));
                }
            } else {
                return false;
            }

            if (auto [ok, arr] = getStringArray(database, "module_path_prefix", s); ok) {
                // Prefix path
                if (arr) {
                    Simulator::prependModulePath(std::move(arrayToStringVector(arr)));
                }
            } else {
                return false;
            }
            if (auto [ok, arr] = getStringArray(database, "module_path_suffix", s); ok) {
                // Suffix path
                if (arr) {
                    Simulator::appendModulePath(std::move(arrayToStringVector(arr)));
                }
            } else {
                return false;
            }
        }

        if (auto database = config["Binaries"].as_table()) {
            if (auto [ok, found, str] = getString(database, "openvaf", s); ok) {
                if (found) {
                    Platform::setOpenVaf(str);
                }
            } else {
                return false;
            }
            
            if (auto [ok, arr] = getStringArray(database, "openvaf_args", s); ok) {
                // Prefix path
                if (arr) {
                    Platform::setOpenVafArgs(std::move(arrayToStringVector(arr)));
                }
            } else {
                return false;
            }

            if (auto [ok, found, str] = getString(database, "python", s); ok) {
                if (found) {
                    Platform::setPythonExecutable(str);
                }
            } else {
                return false;
            }
        }
    } catch (const toml::parse_error& err) {
        std::ostringstream oss;
        oss << "Parsing of file " << filename << " failed:\n"
            << "  " << err.description() << "\n"
            << "  Line: " << err.source().begin.line 
            << ", Column: " << err.source().begin.column;
        
        s.set(Status::Syntax, oss.str());
        return false;
    }
    
    return true;
}

}
