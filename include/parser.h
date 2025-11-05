#ifndef __PARSERDRIVER_DEFINED
#define __PARSERDRIVER_DEFINED

#include <string>
#include <cstddef>
#include <istream>
#include <unordered_map>
#include "value.h"
#include "status.h"
#include "location.h"
#include "parseroutput.h"
#include "filestack.h"
#include "rpneval.h"
#include "common.h"


namespace NAMESPACE {

class Parser {
public:
    Parser(ParserTables& tab) : tab_(tab) {};
    virtual ~Parser();

    Parser           (const Parser&)  = delete;
    Parser           (      Parser&&) = default;
    Parser& operator=(const Parser&)  = delete;
    Parser& operator=(      Parser&&) = default;
    
    // Do not throw, return false on error
    bool parseNetlistFile(FileStackFileIndex fileIndex, Status& s=Status::ignore);
    
    bool parseNetlistString(const std::string& input, Status& s=Status::ignore);
    bool parseNetlistString(const std::string&& input, Status& s=Status::ignore);
    
    // Throw on error
    Rpn parseExpression(const std::string& input, Status& s=Status::ignore);
    Rpn parseExpression(const std::string&& input, Status& s=Status::ignore);

    bool parseParameters(const std::string& input, Status& s=Status::ignore);
    bool parseParameters(const std::string&& input, Status& s=Status::ignore);

private:
    bool netlistParseHelper(std::istream& stream, Status& s=Status::ignore);

    bool exprParseHelper(std::istream& stream, Status& s=Status::ignore);

    ParserTables& tab_;
    Rpn parsedExpression;
    PTParameters parsedParameters;
};

}

#endif
