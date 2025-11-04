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
#include "dflscanner.h"
#include "dflparser.h"
#include "filestack.h"
#include "rpneval.h"
#include "common.h"


namespace NAMESPACE {

class Parser {
public:
    Parser();
    virtual ~Parser();

    Parser           (const Parser&)  = delete;
    Parser           (      Parser&&) = default;
    Parser& operator=(const Parser&)  = delete;
    Parser& operator=(      Parser&&) = default;
    
    bool parseNetlistFile(FileStackFileIndex fileIndex, ParserTables& tab, ParserExtras& extras, Status& s=Status::ignore);
    
    bool parseNetlistString(const std::string& input, ParserTables& tab, ParserExtras& extras, Status& s=Status::ignore);
    bool parseNetlistString(const std::string&& input, ParserTables& tab, ParserExtras& extras, Status& s=Status::ignore);
    
    bool parseExpression(const std::string& input, Status& s=Status::ignore);
    bool parseExpression(const std::string&& input, Status& s=Status::ignore);

    bool parseParameters(const std::string& input, Status& s=Status::ignore);
    bool parseParameters(const std::string&& input, Status& s=Status::ignore);


private:
    bool netlistParseHelper(std::istream& stream, ParserTables& tab, ParserExtras& extras, Status& s=Status::ignore);
};

}

#endif
