#include <cctype>
#include <fstream>
#include <cassert>
#include <string>
#include <sstream>
#include "parser.h"
#include "common.h"


namespace NAMESPACE {

Parser::Parser() {
}
    
Parser::~Parser() {
}

bool Parser::parseNetlistFile(const char* const filename, ParserTables& tab, ParserExtras& extras, Status& s) {
    assert( filename != nullptr );
    
    if (Simulator::fileDebug()) {
        Simulator::dbg() << "Opening file '" << filename << "'.\n";
    }
                    
    auto stackPosition = tab.fileStack().addFile(filename);
    if (stackPosition==FileStack::badFileId) {
        s.set(Status::NotFound, std::string("File '")+filename+"' not found.");
        return false;
    }
    
    std::ifstream in_file(tab.fileStack().canonicalName(stackPosition));
    if(!in_file.good()) {
        s.set(Status::NotFound, std::string("Failed to open file '")+filename+"'.");
        return false;
    }
    return netlistParseHelper(in_file, tab, extras, s);
}

bool Parser::parseNetlistString(const std::string& input, ParserTables& tab, ParserExtras& extras, Status& s) {
    std::istringstream stream;
    stream.str(input);
    tab.fileStack().addStringFile(input);

    return netlistParseHelper(stream, tab, extras, s); 
}

bool Parser::parseNetlistString(const std::string&& input, ParserTables& tab, ParserExtras& extras, Status& s) {
    std::istringstream stream;
    stream.str(std::move(input));
    tab.fileStack().addStringFile(input);

    return netlistParseHelper(stream, tab, extras, s); 
}

bool Parser::netlistParseHelper(std::istream &stream, ParserTables& tab, ParserExtras& extras, Status& s) {
    RpnEvaluator evaluator;
    dflparse::Scanner scanner(&stream, tab, dflparse::Scanner::InputNetlist, s);
    dflparse::Parser parser(scanner, tab, extras, evaluator, s);
    
    const int accept(0);
    if (parser.parse() != accept){
        return false;
    }

    return true;
}

}
