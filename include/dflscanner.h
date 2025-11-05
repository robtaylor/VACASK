#ifndef __DFLSCANNER_DEFINED
#define __DFLSCANNER_DEFINED

#if ! defined(yyFlexLexerOnce)
#include <FlexLexer.h>
#endif

#include <fstream>
#include <sstream>
#include <string>
#include "dflparser.h"
#include "location.h"
#include "status.h"
#include "simulator.h"
#include "common.h"


namespace NAMESPACE::dflparse {

class Scanner : public yyFlexLexer{
public:   
    typedef enum InputType { InputNetlist, InputExpression, InputParameters } InputType; 

    Scanner(std::istream* in, ParserTables& tab, InputType inputType, FileStackFileIndex fileIndex, Status& s=Status::ignore) 
        : yyFlexLexer(in), tables(tab), inputType(inputType), atBeginning(true), fileStackPosition (fileIndex), status_(s), 
          inParen(0), inBracket(0), inControl(false) {};
    virtual ~Scanner() { cleanup(); };

    Scanner           (const Scanner&)  = delete;
    Scanner           (      Scanner&&) = default;
    Scanner& operator=(const Scanner&)  = delete;
    Scanner& operator=(      Scanner&&) = default;

    // Cleanup 
    void cleanup() { streamStack.clear(); locationStack.clear(); }; 
    
    // Get rid of override virtual function warning
    using FlexLexer::yylex;

    virtual int yylex(NAMESPACE::dflparse::Parser::semantic_type * const lval, NAMESPACE::dflparse::Parser::location_type *location );
    
    // YY_DECL defined in dfllexer.l
    // Method body created by flex in dfllexer.cpp

    std::ifstream* pushStream(const std::string& fileName, Location& loc) {
        streamStack.emplace_back();
        auto& s = streamStack.back();
        s.open(fileName);
        locationStack.push_back(loc);
        if (s.good()) {
            return &(streamStack.back());
        } else {
            return nullptr;
        }
    };

    const char* skipLeadingWhitespace(const char* txt) {
        const char* ch;
        for(ch=txt; *ch; ch++) {
            if (*ch==' ' || *ch=='\t' || *ch=='\n' || *ch=='\r') {
                continue;
            }
            break;
        }
        return ch;
    };

    Location popStream() { streamStack.pop_back(); Location l=locationStack.back(); locationStack.pop_back(); return l; };
    
    void error( const Parser::location_type &l, const std::string &err_message ) {
        status_.set(Status::Syntax, err_message);
        status_.extend(l.loc());
    };

private:
    // yyval ptr
    Parser::semantic_type *yylval = nullptr;
    ParserTables& tables;
    FileStackFileIndex fileStackPosition;
    std::string sbuf;
    std::string marker;
    Position stringStart;
    std::vector<std::ifstream> streamStack;
    std::vector<Location> locationStack;
    Status& status_;
    std::string section;
    size_t inParen;
    size_t inBracket;
    bool inControl;
    InputType inputType;
    bool atBeginning;
};

}

#endif
