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

bool Parser::parseNetlistFile(FileStackFileIndex stackPosition, ParserTables& tab, ParserExtras& extras, Status& s) {
    auto t0 = Accounting::wclk();
    tab.accounting().acctNew.parse++;
    
    auto& fileName = tab.fileStack().canonicalName(stackPosition);

    if (Simulator::fileDebug()) {
        Simulator::dbg() << "Reading file '" << fileName << "'.\n";
    }
                    
    std::ifstream in_file(tab.fileStack().canonicalName(stackPosition));
    if(!in_file.good()) {
        s.set(Status::NotFound, std::string("Failed to open file '")+fileName+"'.");
        tab.accounting().acctNew.tparse += Accounting::wclkDelta(t0);
        return false;
    }
    auto st = netlistParseHelper(in_file, tab, extras, s);
    tab.accounting().acctNew.tparse += Accounting::wclkDelta(t0);
    return st;
}

bool Parser::parseNetlistString(const std::string& input, ParserTables& tab, ParserExtras& extras, Status& s) {
    auto t0 = Accounting::wclk();
    tab.accounting().acctNew.parse++;
    
    std::istringstream stream;
    stream.str(input);
    tab.fileStack().addStringFile(input);

    auto st = netlistParseHelper(stream, tab, extras, s); 
    tab.accounting().acctNew.tparse += Accounting::wclkDelta(t0);
    return st;
}

bool Parser::parseNetlistString(const std::string&& input, ParserTables& tab, ParserExtras& extras, Status& s) {
    auto t0 = Accounting::wclk();
    std::istringstream stream;
    stream.str(std::move(input));
    tab.fileStack().addStringFile(input);

    auto st = netlistParseHelper(stream, tab, extras, s); 
    tab.accounting().acctNew.tparse += Accounting::wclkDelta(t0);
    return st;
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
