#include <cctype>
#include <fstream>
#include <cassert>
#include <string>
#include <sstream>
#include "parser.h"
#include "dflscanner.h"
#include "dflparser.h"
#include "common.h"


namespace NAMESPACE {

Parser::~Parser() {
}

bool Parser::parseNetlistFile(FileStackFileIndex stackPosition, Status& s) {
    auto t0 = Accounting::wclk();
    tab_.accounting().acctNew.parse++;
    
    auto& fileName = tab_.fileStack().canonicalName(stackPosition);

    if (Simulator::fileDebug()) {
        Simulator::dbg() << "Reading file '" << fileName << "'.\n";
    }
                    
    std::ifstream in_file(tab_.fileStack().canonicalName(stackPosition));
    if(!in_file.good()) {
        s.set(Status::NotFound, std::string("Failed to open file '")+fileName+"'.");
        tab_.accounting().acctNew.tparse += Accounting::wclkDelta(t0);
        return false;
    }
    auto st = netlistParseHelper(in_file, stackPosition, s);
    tab_.accounting().acctNew.tparse += Accounting::wclkDelta(t0);
    return st;
}

bool Parser::parseNetlistString(const std::string& input, Status& s) {
    auto t0 = Accounting::wclk();
    tab_.accounting().acctNew.parse++;
    
    // Add to filestack (makes a copy)
    auto pos = tab_.fileStack().addStringFile(input);

    // Create a stream, copy string into stream
    std::istringstream stream;
    stream.str(input);
    
    auto st = netlistParseHelper(stream, pos, s); 
    tab_.accounting().acctNew.tparse += Accounting::wclkDelta(t0);
    return st;
}

bool Parser::parseNetlistString(const std::string&& input, Status& s) {
    auto t0 = Accounting::wclk();
    
    // Add to filestack (makes a copy)
    auto pos = tab_.fileStack().addStringFile(input);

    // Create a stream, move string into stream
    std::istringstream stream;
    stream.str(std::move(input));
    
    auto st = netlistParseHelper(stream, pos, s); 
    tab_.accounting().acctNew.tparse += Accounting::wclkDelta(t0);
    return st;
}

bool Parser::netlistParseHelper(std::istream &stream, FileStackFileIndex pos, Status& s) {
    RpnEvaluator evaluator;

    dflparse::Scanner scanner(&stream, tab_, dflparse::Scanner::InputNetlist, pos, s);
    dflparse::Parser parser(scanner, tab_, nullptr, nullptr, evaluator, s);
    
    const int accept(0);
    if (parser.parse() != accept){
        return false;
    }

    return true;
}

bool Parser::exprParseHelper(std::istream &stream, FileStackFileIndex pos, Status& s) {
    RpnEvaluator evaluator;

    dflparse::Scanner scanner(&stream, tab_, dflparse::Scanner::InputExpression, pos, s);
    dflparse::Parser parser(scanner, tab_, &parsedExpression, nullptr, evaluator, s);
    
    const int accept(0);
    if (parser.parse() != accept){
        return false;
    }

    return true;
}

Rpn Parser::parseExpression(const std::string& input, Status& s) {
    auto t0 = Accounting::wclk();

    // Add to filestack (makes a copy)
    auto pos = tab_.fileStack().addStringFile(input);

    // Create a stream, copy string into stream
    std::istringstream stream;
    stream.str(input);

    Status tmps;
    auto st = exprParseHelper(stream, pos, tmps); 
    
    tab_.accounting().acctNew.tparse += Accounting::wclkDelta(t0);

    // Throw on failure
    if (!st) {
        s.set(tmps);
        Simulator::err() << tmps.message();
        std::runtime_error(std::string("Failed to parse expression '")+stream.str()+"'"); 
    }

    return std::move(parsedExpression);
}

Rpn Parser::parseExpression(const std::string&& input, Status& s) {
    auto t0 = Accounting::wclk();

    // Add to filestack (makes a copy)
    auto pos = tab_.fileStack().addStringFile(input);

    // Create a stream, move string into stream
    std::istringstream stream;
    stream.str(std::move(input));

    Status tmps;
    auto st = exprParseHelper(stream, pos, tmps); 
    
    tab_.accounting().acctNew.tparse += Accounting::wclkDelta(t0);

    // Throw on failure
    // TODO: collect status even if ignore is passed
    //       forward only if not ignore
    if (!st) {
        s.set(tmps);
        Simulator::err() << tmps.message();
        throw std::runtime_error(std::string("Failed to parse expression '")+stream.str()+"'"); 
    }

    return std::move(parsedExpression);
}

bool Parser::parametersParseHelper(std::istream &stream, FileStackFileIndex pos, Status& s) {
    RpnEvaluator evaluator;

    dflparse::Scanner scanner(&stream, tab_, dflparse::Scanner::InputParameters, pos, s);
    dflparse::Parser parser(scanner, tab_, nullptr, &parsedParameters, evaluator, s);
    
    const int accept(0);
    if (parser.parse() != accept){
        return false;
    }

    return true;
}

PTParameters Parser::parseParameters(const std::string& input, Status& s) {
    auto t0 = Accounting::wclk();

    // Add to filestack (makes a copy)
    auto pos = tab_.fileStack().addStringFile(input);

    // Create a stream, copy string into stream
    std::istringstream stream;
    stream.str(input);

    Status tmps;
    auto st = parametersParseHelper(stream, pos, tmps); 
    
    tab_.accounting().acctNew.tparse += Accounting::wclkDelta(t0);

    // Throw on failure
    if (!st) {
        s.set(tmps);
        Simulator::err() << tmps.message();
        std::runtime_error(std::string("Failed to parse parameters '")+stream.str()+"'"); 
    }

    return std::move(parsedParameters);
}

PTParameters Parser::parseParameters(const std::string&& input, Status& s) {
    auto t0 = Accounting::wclk();

    // Add to filestack (makes a copy)
    auto pos = tab_.fileStack().addStringFile(input);

    // Create a stream, move string into stream
    std::istringstream stream;
    stream.str(std::move(input));

    Status tmps;
    auto st = parametersParseHelper(stream, pos, tmps); 
    
    tab_.accounting().acctNew.tparse += Accounting::wclkDelta(t0);

    // Throw on failure
    if (!st) {
        s.set(tmps);
        Simulator::err() << tmps.message();
        throw std::runtime_error(std::string("Failed to parse expression '")+stream.str()+"'"); 
    }

    return std::move(parsedParameters);
}

}
