%skeleton "lalr1.cc"

// Minimum version of Bison
%require  "3.3"

// Extra output file with description of states
// %verbose 

// Turn on parser instrumentation for tracing
// %define parse.trace

// Verbose parse errors
%define parse.error verbose

// Write a header file with token definitions (will change to %header)
%defines

// Namespace to use for the parser
%define api.namespace {NAMESPACE::dflparse}

// Name of the parser class
%define api.parser.class {Parser}

// Code required for the value and location types
// Goes to the top of the parser include file
%code requires{
// #define YYDEBUG 1
#include "value.h"
#include "parseroutput.h"
#include "rpneval.h"
#include "rpnexpr.h"
#include "location.h"
#include <utility>
#include <tuple>
#include <iostream>
#include "common.h"

using namespace NAMESPACE;

namespace NAMESPACE::dflparse {
    class Scanner;
}

namespace NAMESPACE {
    // Stuff from library namespace
    static Id saveCmd = Id::createStatic("save");
}

// The following definitions is missing when %locations isn't used
# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

struct pexpr {
    Id id; 
    Rpn expr; 
    Location loc;
};

struct paramlist {
    PTParameters params;
    std::unordered_map<Id,Location> locations;
};

struct sweeps {
    std::vector<PTSweep> sweeps;
    std::unordered_map<Id,Location> locations;
};

typedef struct subckt {
    bool isToplevel {true};
    bool hasControlBlock {false};
    PTSubcircuitDefinition def;
    PTParameters parameters;
    std::unordered_map<Id,Location> paramLoc;
} subckt;

}

// Additional arguments to parser class constructor
%parse-param {NAMESPACE::dflparse::Scanner& scanner}
// %parse-param {ParserDriver& driver}
%parse-param {ParserTables& tables}
%parse-param {Rpn* expressionPtr}
%parse-param {PTParameters* parametersPtr}
%parse-param {RpnEvaluator& evaluator}
%parse-param {Status &status}

// Add this code to the beginning of parser implementation
%code{
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <fstream>

#include "dflscanner.h"

#undef yylex 
#define yylex scanner.yylex
}
 
// Type of semantic value, similar to union, allows any c++ object type
%define api.value.type variant

// Include runtime assertions to check for invalid use
// In c++ parsers it uses runtime type information (RTTI)
%define parse.assert

// Generate code for location tracking
%locations

// Location type
%define api.location.type {Location}

// Tokens, semantic types, and token names for error reporting
%token               END    0     "end of file"
%token <std::string> TITLE        "circuit title"

%token               MODEL        "model"
%token               GLOBAL       "global"
%token               GROUND       "ground"
%token               OPTIONS      "options"
%token               LOAD         "load"
%token               SUBCKT       "subckt"
%token               ENDS         "ends"
%token               PARAMETERS   "parameters"
%token               EMBED        "embed"
%token               CONTROL      "control"
%token               ENDC         "endc"
%token               SAVE         "save"
%token               SWEEP        "sweep"
%token               ANALYSIS     "analysis"

%token <Id>          IDENTIFIER   "identifier"

%token <Int>         INTEGER      "integer"
%token <Int>         HEXINTEGER   "hexadecimal integer"
%token <Real>        FLOAT        "floating point number"
%token <std::string> STRING       "string literal"

%token               PLUS         "+"
%token               MINUS        "-"
%token               TIMES        "*"
%token               DIVIDE       "/"
%token               POWER        "**"
%token               LBRACKET     "["
%token               RBRACKET     "]"
%token               LPAREN       "("
%token               RPAREN       ")"
%token               ASSIGN       "="
%token               COMMA        ","
%token               GREATER      ">"
%token               LESS         "<"
%token               GREATEREQ    ">="
%token               LESSEQ       "<="
%token               EQUAL        "=="
%token               NOTEQUAL     "!="
%token               QUESTION     "?"
%token               AND          "&&"
%token               OR           "||"
%token               NOT          "!"
%token               BITAND       "&"
%token               BITOR        "|"
%token               BITNOT       "~"
%token               BITEXOR      "^"
%token               BITSHIFTR    ">>"
%token               BITSHIFTL    "<<"

%token               COLON        ":"
%token               SEMICOLON    ";"
%token               RIGHTARROW   "->"
// Not allowed due to conflict with <-5 (could be <- 5 or < -5)
// %token               LEFTARROW    "<-" 

%token               NEWLINE      "newline"

%token               INNETLIST    "input netlist"
%token               INEXPR       "input expression"
%token               INPARAMS     "input parameters"

%token               BLKIF        "@if"
%token               BLKELSEIF    "@elseif"
%token               BLKELSE      "@else"
%token               BLKENDIF     "@endif"


// Operator associativity and precedence, lowest first
%right QUESTION
%left OR
%left AND
%left BITOR
%left BITEXOR
%left BITAND
%left EQUAL NOTEQUAL 
%left LESS GREATER GREATEREQ LESSEQ
%left BITSHIFTL BITSHIFTR
%left PLUS MINUS
%left TIMES DIVIDE
%right POWER
%right NEG NOT BITNOT
%left LPAREN RPAREN LBRACKET RBRACKET

// exprlist e  e,e
// bracketlist1 [ [1 
// commabracketlist1 [, [1, [1,2,3
// semicolonbracketlist1 [; [1; [1;2;3

// Nonterminal symbols
%type <Int>                             intnum
%type <Id>                              terminal
%type <PTIdentifierList>                terminal_list global ground keywords
%type <std::vector<Rpn>>                exprlist semexprlist colexprlist
%type <Value>                           value
%type <struct pexpr>                    parameter_expression
%type <struct paramlist>                parameter_list opt_broken_parameter_list subcktparameters 
%type <PTInstance>                      instance
%type <PTModel>                         model
%type <PTSubcircuitDefinition>          subckt
%type <struct subckt>                   subckt_build
%type <PTBlockSequence>                 condblock_build condblock
%type <PTLoad>                          load
%type <Rpn>                             expr
%type <Id>                              savestr
%type <std::vector<Id>>                 savestrlist
%type <PTSave>                          savecmd
%type <std::vector<PTSave>>             savecmd_list saves
%type <PTEmbed>                         embed
%type <Int>                             control_block_build control_block
%type <PTCommand>                       command
%type <struct sweeps>                   sweeps
%type <PTAnalysis>                      analysis pre_analysis analysis_with_params 


// Rules
%%

output
  : INNETLIST subckt_build END {
    // Toplevel circuit definition
    $2.def.add(std::move($2.parameters));
    tables.setDefaultSubDef(std::move($2.def));
    tables.defaultGround();
    // Verify tables (basic level 0 verifications)
    if (!(tables.verify(status))) {
        YYERROR;
    }
  }
  | INEXPR expr END {
    // Parse an expression
    *expressionPtr = std::move($2);
  }
  | INPARAMS opt_broken_parameter_list END {
    // Parse a parameters list
    *parametersPtr = std::move($2.params);
  }

// Toplevel netlist and subcircuit definition
subckt_build
  : TITLE NEWLINE {
    // Toplevel netlist start
    tables.setTitle(std::move($1));
    $$.def = std::move(PTSubcircuitDefinition(
        Id(), // By default the name (Id) of the toplevel definition is not valid (empty)
        std::move(PTIdentifierList()), 
        @1.loc()
    ));
  }
  | SUBCKT IDENTIFIER LPAREN RPAREN NEWLINE {
    // Subcircuit definition start, no terminals
    $$.def = std::move(PTSubcircuitDefinition(
        $2, 
        PTIdentifierList(), 
        @1.loc()
    ));
    // This is not the toplevel definition
    $$.isToplevel = false;
  }
  | SUBCKT IDENTIFIER LPAREN terminal_list RPAREN NEWLINE {
    // Subcircuit definition start, with terminals
    $$.def = std::move(PTSubcircuitDefinition(
        $2, 
        std::move($4), 
        @1.loc()
    ));
    // This is not the toplevel definition
    $$.isToplevel = false;
  }
  | subckt_build NEWLINE {
    // Skip NEWLINE
    $$ = std::move($1);
  }
  | subckt_build subcktparameters {
    // parameters
    for(auto it=$2.locations.begin(); it!=$2.locations.end(); ++it) {
        auto fdit = $1.paramLoc.find(it->first);
        if (fdit!=$1.paramLoc.end()) {
            status.set(Status::Redefinition, "Parameter redefinition.");
            status.extend(it->second.loc());
            status.extend("Parameter first defined here.");
            status.extend(fdit->second.loc());
            YYERROR;
        }
    }
    $$ = std::move($1);
    $$.parameters.add(std::move($2.params));
    $$.paramLoc.merge(std::move($2.locations));
  }
  | subckt_build model {
    $$ = std::move($1);
    $$.def.add(std::move($2));
  }
  | subckt_build instance {
    $$ = std::move($1);
    $$.def.add(std::move($2));
  }
  | subckt_build condblock {
    $$ = std::move($1);
    $$.def.add(std::move($2));
  }
  | subckt_build subckt {
    // subcircuit definition, not allowed inside other subcircuit definitions
    if (!$1.isToplevel) {
        status.set(Status::Syntax, "Nested subcircuit definitions are not allowed.");
        status.extend(@2.loc());
        YYERROR;
    }
    $$ = std::move($1);
    $$.def.add(std::move($2));
  }
  | subckt_build global {
    // global nodes, not allowed in subcircuit definition
    if (!$$.isToplevel) {
        status.set(Status::Syntax, "Global nodes allowed only in toplevel circuit.");
        status.extend(@2.loc());
        YYERROR;
    }
    $$ = std::move($1);
    for(auto it=$2.begin(); it!=$2.end(); ++it) {
        tables.addGlobal(std::move(*it));
    }
  }
  | subckt_build ground {
    // ground nodes, not allowed in subcircuit definition
    if (!$$.isToplevel) {
        status.set(Status::Syntax, "Ground nodes allowed only in toplevel circuit.");
        status.extend(@2.loc());
        YYERROR;
    }
    $$ = std::move($1);
    for(auto it=$2.begin(); it!=$2.end(); ++it) {
        tables.addGround(std::move(*it));
    }
  }
  | subckt_build load {
    // load, not allowed in subcircuit definition
    if (!$$.isToplevel) {
        status.set(Status::Syntax, "Load allowed only in toplevel circuit.");
        status.extend(@2.loc());
        YYERROR;
    }
    $$ = std::move($1);
    tables.add(std::move($2));
  }
  | subckt_build embed {
    // embed, not allowed in subcircuit definition
    if (!$$.isToplevel) {
        status.set(Status::Syntax, "Embed allowed only in toplevel circuit.");
        status.extend(@2.loc());
        YYERROR;
    }
    $$ = std::move($1);
    tables.add(std::move($2));
  }
  | subckt_build control_block {
    // control block, allow this for toplevel definition only
    if (!$1.isToplevel) {
        status.set(Status::Syntax, "Control block is not allowed inside subcircuit definition.");
        status.extend(@2.loc());
        YYERROR;
    }
    // Only one control block is allowed
    if ($1.hasControlBlock) {
        status.set(Status::Syntax, "Only one control block is allowed.");
        status.extend(@2.loc());
        YYERROR;
    }
    $$ = std::move($1);
    $$.hasControlBlock = true;
  }
  
subckt
  : subckt_build ENDS NEWLINE {
    // Subcircuit definition
    $1.def.add(std::move($1.parameters)); 
    $$ = std::move($1.def);
  }

// Conditional block building
condblock_build 
  : BLKIF expr NEWLINE {
    PTBlock blk;
    $$.add(std::move($2), std::move(blk), @1.loc());
  }
  | condblock_build instance {
    $$ = std::move($1);
    $$.back().add(std::move($2));
  }
  | condblock_build model {
    $$ = std::move($1);
    $$.back().add(std::move($2));
  }
  | condblock_build condblock {
    $$ = std::move($1);
    $$.back().add(std::move($2));
  }
  | condblock_build BLKELSEIF expr NEWLINE {
    $$ = std::move($1);
    PTBlock blk;
    $$.add(std::move($3), std::move(blk), @2.loc());
  }
  | condblock_build BLKELSE NEWLINE {
    $$ = std::move($1);
    PTBlock blk;
    $$.add(std::move(Rpn()), std::move(blk), @2.loc());
  }
  | condblock_build NEWLINE {
    $$ = std::move($1);
  }
;

// End of conditional block
condblock
  : condblock_build BLKENDIF NEWLINE {
    $$ = std::move($1);
  }
;

terminal
  : IDENTIFIER { 
    $$ = $1;
  }
  | INTEGER { 
    $$ = Id(std::to_string($1)); 
  }

terminal_list
  : terminal { 
    $$.push_back(PTParsedIdentifier(std::move($1), @1.loc())); 
  }
  | terminal_list terminal { 
    $$ = std::move($1); 
    $$.push_back(PTParsedIdentifier(std::move($2), @2.loc())); 
  } 

intnum
  : INTEGER { $$ = $1; } 
  | HEXINTEGER { $$ = $1; } 

value
  : intnum { $$ = Int($1); }
  | FLOAT { $$ = Real($1); }
  | STRING { $$ = Value($1); }

expr
  : value { $$.extend(std::move($1), @1.loc()); }
  | IDENTIFIER { 
    $$.extend(Rpn::Identifier(std::move($1)), @1.loc()); 
  }
  | expr PLUS expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpPlus), @2.loc()); 
  }
  | expr MINUS expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpMinus), @2.loc()); 
  }
  | expr TIMES expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpTimes), @2.loc());  
  }
  | expr DIVIDE expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpDivide), @2.loc());  
  }
  | expr POWER expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpPower), @2.loc());  
  }
  | expr EQUAL expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpEqual), @2.loc());  
  }
  | expr NOTEQUAL expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpEqual), @2.loc());  
  }
  | expr LESS expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpLess), @2.loc());  
  }
  | expr LESSEQ expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpLessEq), @2.loc());  
  }
  | expr GREATER expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpGreater), @2.loc());  
  }
  | expr GREATEREQ expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpGreaterEq), @2.loc());  
  }
  | expr BITAND expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpBitAnd), @2.loc());  
  }
  | expr BITOR expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpBitOr), @2.loc());  
  }
  | expr BITEXOR expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpBitExor), @2.loc());  
  }
  | expr BITSHIFTR expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpBitShiftR), @2.loc());  
  }
  | expr BITSHIFTL expr { 
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpBitShiftL), @2.loc());  
  }
  | expr AND expr { 
    // short circuit (a && b) translation to RPN
    //         a
    //         makeboolean
    //         branchiffalse end
    //         b
    //         makeboolean
    //         op(and) // does nothing during execution, needed by formatting
    //   end:
    auto tailLen = $3.size();
    auto needsConversion = !$3.endsWithMakeBoolean();
    $$.extend(std::move($1)); 
    $$.extend(Rpn::MakeBoolean(), @2.loc());
    $$.extend(Rpn::Branch($3.size()+(needsConversion?1:0)+2, Rpn::BrFalse|Rpn::BrKeepOnBranch|Rpn::BrHidden), @2.loc());
    $$.extend(std::move($3)); 
    if (needsConversion) {
      $$.extend(Rpn::MakeBoolean(), @2.loc());
    }
    $$.extend(Rpn::Op(Rpn::OpAnd), @2.loc());  
  }
  | expr OR expr { 
    // short circuit (a &|| b) translation to RPN
    //         a
    //         makeboolean
    //         branchiftrue end
    //         b
    //         makeboolean
    //         op(or) // does nothing during execution, needed by formatting
    //   end:
    auto tailLen = $3.size();
    auto needsConversion = !$3.endsWithMakeBoolean();
    $$.extend(std::move($1)); 
    $$.extend(Rpn::MakeBoolean(), @2.loc());
    $$.extend(Rpn::Branch($3.size()+(needsConversion?1:0)+2, Rpn::BrKeepOnBranch|Rpn::BrHidden), @2.loc());
    $$.extend(std::move($3)); 
    if (needsConversion) {
      $$.extend(Rpn::MakeBoolean(), @2.loc());
    }
    $$.extend(Rpn::Op(Rpn::OpOr), @2.loc());
  }
  | expr QUESTION expr COLON expr %prec QUESTION {
    // Ternary operator a?b:c, translation to RPM
    //        a
    //        branchiffalse false // +3
    //        b
    //        jump end // +2
    // false: c
    // end:   op(question) // does nothing during execution, needed by formatting
    $$.extend(std::move($1)); 
    $$.extend(Rpn::Branch(3, Rpn::BrFalse|Rpn::BrHidden), @2.loc());
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Jump(2, Rpn::BrHidden), @2.loc());
    $$.extend(std::move($5)); 
    $$.extend(Rpn::Op(Rpn::OpQuestion), @2.loc());
  }
  | BITNOT expr { 
    $$.extend(std::move($2)); 
    $$.extend(Rpn::Op(Rpn::OpBitNot), @1.loc());  
  }
  | NOT expr { 
    $$.extend(std::move($2)); 
    $$.extend(Rpn::Op(Rpn::OpNot), @1.loc());  
  }
  | MINUS expr %prec NEG { 
    $$.extend(std::move($2)); 
    $$.extend(Rpn::Op(Rpn::OpUMinus), @1.loc());  
  }
  | PLUS expr %prec NEG { 
    $$.extend(std::move($2)); 
  }
  | LPAREN expr RPAREN{ 
    $$.extend(std::move($2)); 
  }
  | IDENTIFIER LPAREN RPAREN { 
    // Function call, no arguments
    $$.extend(Rpn::FunctionCall(std::move($1), 0), @1.loc()); 
  }
  | IDENTIFIER LPAREN exprlist RPAREN { 
    // Function call with arguments
    for(Rpn::Arity i=0; i<$3.size(); i++) {
        $$.extend(std::move($3[i])); 
    }
    $$.extend(Rpn::FunctionCall(std::move($1), $3.size()), @1.loc()); 
  }
  | LBRACKET RBRACKET {
    // Empty vector of type Int
    $$.extend(Rpn::PackVec(0), @1.loc()); 
  }
  | LBRACKET COMMA RBRACKET {
    // Empty vector of type Int
    $$.extend(Rpn::PackVec(0), @1.loc()); 
  }
  | LBRACKET exprlist RBRACKET {
    // Pack values in a vector, flatten lists into a vector
    // [,] is an empty Int vector
    for(Rpn::Arity i=0; i<$2.size(); i++) {
        $$.extend(std::move($2[i])); 
    }
    $$.extend(Rpn::PackVec($2.size()), @1.loc()); 
  }
  | LBRACKET SEMICOLON RBRACKET {
    // Empty list
    $$.extend(Rpn::PackList(0), @1.loc()); 
  }
  | LBRACKET expr SEMICOLON RBRACKET {
    // List with single element
    // Pack values in a list, keep members that are lists themselves intact
    // This produces a list of lists of ...
    $$.extend(std::move($2)); 
    $$.extend(Rpn::PackList(1), @1.loc()); 
  }
  | LBRACKET semexprlist RBRACKET {
    // List with two or more elements
    // Pack values in a list, keep members that are lists themselves intact
    // This produces a list of lists of ...
    for(Rpn::Arity i=0; i<$2.size(); i++) {
        $$.extend(std::move($2[i])); 
    }
    $$.extend(Rpn::PackList($2.size()), @1.loc()); 
  }
  | LBRACKET COLON RBRACKET {
    // Empty list
    $$.extend(Rpn::PackList(0), @1.loc()); 
  }
  | LBRACKET expr COLON RBRACKET {
    // List with single element upacked
    // Merge scalars and lists in one list
    $$.extend(std::move($2)); 
    $$.extend(Rpn::MergeList(1), @1.loc()); 
  }
  | LBRACKET colexprlist RBRACKET {
    // List with two or more elements unpacked
    // Merge scalars and lists in one list
    for(Rpn::Arity i=0; i<$2.size(); i++) {
        $$.extend(std::move($2[i])); 
    }
    $$.extend(Rpn::MergeList($2.size()), @1.loc()); 
  }
  
  | expr LBRACKET expr RBRACKET {
    // Vector and list selector
    $$.extend(std::move($1)); 
    $$.extend(std::move($3)); 
    $$.extend(Rpn::Op(Rpn::OpSelect), @2.loc()); 
  }

// Comma separated exprlist is always in parentheses or brackets, 
// no need to handle NEWLINE. 
// List has always at least one expression. 
exprlist
  : expr {
    // Single expresion
    $$.push_back(std::move($1));
  }
  | exprlist COMMA expr {
    // Multiple expressions
    $$ = std::move($1);
    $$.push_back(std::move($3));
  }

// Semicolon separated expression list is always in brackets, 
// no need to handle NEWLINE. 
// List has always at least two expressions. 
semexprlist
  : expr SEMICOLON expr {
    $$.push_back(std::move($1));
    $$.push_back(std::move($3));
  }
  | semexprlist SEMICOLON expr {
    $$ = std::move($1);
    $$.push_back(std::move($3));
  }

// Colon separated expression list is always in brackets, 
// no need to handle NEWLINE. 
// List has always at least two expressions. 
colexprlist
  : expr COLON expr {
    $$.push_back(std::move($1));
    $$.push_back(std::move($3));
  }
  | colexprlist COLON expr {
    $$ = std::move($1);
    $$.push_back(std::move($3));
  }

parameter_expression
  : IDENTIFIER ASSIGN expr { 
    $$.id = $1;
    $$.expr = std::move($3);
    $$.loc = @1;
  }

parameter_list
  : parameter_expression {
    $$.locations[$1.id] = $1.loc;
    if (evaluator.isConstant($1.expr)) {
        Value v;
        if (!evaluator.evaluate($1.expr, v, status)) {
            YYERROR;
        }
        $$.params.add(PTParameterValue($1.id, std::move(v), @1.loc()));
        auto dump = std::move($1.expr);
    } else {
        $$.params.add(PTParameterExpression($1.id, std::move($1.expr), @1.loc()));
    }
  }
  | parameter_list parameter_expression {
    $$ = std::move($1);
    auto it = $$.locations.find($2.id);
    if (it!=$$.locations.end()) {
        status.set(Status::Redefinition, "Parameter redefinition.");
        status.extend(@2.loc());
        status.extend("Parameter first defined here.");
        status.extend(it->second.loc());
        YYERROR;
    }
    $$.locations[$2.id] = $2.loc;
    if (evaluator.isConstant($2.expr)) {
        Value v;
        if (!evaluator.evaluate($2.expr, v, status)) {
            YYERROR;
        }
        $$.params.add(PTParameterValue($2.id, std::move(v), @2.loc()));
        auto dump = std::move($2.expr);
    } else {
        $$.params.add(PTParameterExpression($2.id, std::move($2.expr), @2.loc()));
    }
  }

opt_broken_parameter_list
  : parameter_list {
    $$ = std::move($1);
  }
  | LPAREN parameter_list RPAREN {
    $$ = std::move($2);
  }

instance
  : IDENTIFIER LPAREN RPAREN IDENTIFIER NEWLINE {
    // No terminals, no parameters
    $$ = std::move(PTInstance(
        $1, 
        $4, 
        PTIdentifierList(), 
        @1.loc()
    ));
  }
  | IDENTIFIER LPAREN terminal_list RPAREN IDENTIFIER NEWLINE {
    // Terminals, no parameters
    $$ = std::move(PTInstance(
        $1, 
        $5, 
        std::move($3), 
        @1.loc()
    ));
  }
  | IDENTIFIER LPAREN RPAREN IDENTIFIER opt_broken_parameter_list NEWLINE {
    // No terminals, parameters
    $$ = std::move(PTInstance(
        $1, 
        $4, 
        PTIdentifierList(), 
        std::move($5.params), 
        @1.loc()
    ));
  }
  | IDENTIFIER LPAREN terminal_list RPAREN IDENTIFIER opt_broken_parameter_list NEWLINE {
    // Terminals, parameters
    $$ = std::move(PTInstance(
        $1, 
        $5, 
        std::move($3), 
        std::move($6.params), 
        @1.loc()
    ));
  }
  
model
  : MODEL IDENTIFIER IDENTIFIER NEWLINE {
    $$ = std::move(PTModel(
        $2, 
        $3, 
        @1.loc()
    ));
    $$.add(std::move(PTParameters()));
  }
  | MODEL IDENTIFIER IDENTIFIER opt_broken_parameter_list NEWLINE {
    $$ = std::move(PTModel(
        $2, 
        $3, 
        std::move($4.params), 
        @1.loc()
    ));
  }

subcktparameters
  : PARAMETERS opt_broken_parameter_list NEWLINE {
    $$ = std::move($2);
  }

embed
  : EMBED STRING STRING {
    $$ = std::move(PTEmbed(std::move($2), std::move($3), @1.loc())); 
  }

savestr
  : terminal {
    $$ = $1;
  }
  | STRING {
    $$ = $1; 
  }

savestrlist
  : savestr {
    $$.push_back(std::move($1));
  }
  | savestrlist COMMA savestr {
    $$ = std::move($1);
    $$.push_back(std::move($3));
  }

savecmd
  : IDENTIFIER { // LPAREN RPAREN {
    $$ = std::move(PTSave($1, @1.loc()));
  }
  | IDENTIFIER LPAREN savestrlist RPAREN {
    if ($3.size()>2) {
        status.set(Status::BadArguments, "Save directive has too many arguments.");
        status.extend(@1.loc());
        YYERROR;
    } else if ($3.size()==2) {
        $$ = std::move(PTSave($1, $3[0], $3[1], @1.loc()));
    } else {
        $$ = std::move(PTSave($1, $3[0], @1.loc()));
    }
  }

savecmd_list 
  : savecmd {
    $$.push_back(std::move($1));
  }
  | savecmd_list savecmd {
    $$ = std::move($1);
    $$.push_back(std::move($2));
  }

saves
  : savecmd_list {
    $$ = std::move($1);
  }
  | LPAREN savecmd_list RPAREN {
    $$ = std::move($2);
  }

global
  : GLOBAL terminal_list NEWLINE {  
    $$ = std::move($2); 
  }
  | GLOBAL LPAREN terminal_list RPAREN NEWLINE {  
    $$ = std::move($3); 
  }

ground
  : GROUND terminal_list NEWLINE {  
    $$ = std::move($2); 
  }
  | GROUND LPAREN terminal_list RPAREN NEWLINE {  
    $$ = std::move($3); 
  }

load
  : LOAD STRING NEWLINE {
    $$ = std::move(PTLoad($2, @1.loc()));
  }
  | LOAD STRING opt_broken_parameter_list NEWLINE {
    if ($3.params.expressionCount()>0) {
      status.set(Status::BadArguments, "Only constant expressions are allowed here.");
      status.extend($3.params.expressions()[0].location());
      YYERROR;
    }
    $$ = std::move(PTLoad(std::move($2), std::move($3.params), @1.loc()));
  }
  
sweeps
  : SWEEP IDENTIFIER opt_broken_parameter_list {
    Id id = $2;
    $$.sweeps.push_back(PTSweep(id, std::move($3.params), @1.loc()));
    $$.locations.insert({id, @1});
  }
  | sweeps NEWLINE {
    $$ = std::move($1);
  }
  | sweeps SWEEP IDENTIFIER opt_broken_parameter_list {
    $$ = std::move($1);
    Id id = $3;
    auto [it, inserted] = $$.locations.insert({id, @1});
    if (!inserted) {
        status.set(Status::Redefinition, "Sweep does not have a unique name.");
        status.extend(@2.loc());
        status.extend("The name was first used here.");
        status.extend(it->second.loc());
        YYERROR;
    }
    $$.sweeps.push_back(PTSweep(id, std::move($4.params), @2.loc()));
  }

pre_analysis
  : ANALYSIS IDENTIFIER IDENTIFIER {
    $$ = std::move(PTAnalysis($2, $3, @1.loc()));
  }
  | sweeps ANALYSIS IDENTIFIER IDENTIFIER {
    $$ = std::move(PTAnalysis($3, $4, @2.loc()));
    $$.add(std::move($1.sweeps));
  }

analysis_with_params
  : pre_analysis opt_broken_parameter_list {
    $$ = std::move($1);
    $$.add(std::move($2.params));
  }

analysis
  : pre_analysis {
    $$ = std::move($1);
  }
  | analysis_with_params {
    $$ = std::move($1);
  }

keywords
  : IDENTIFIER {
    $$.push_back(PTParsedIdentifier(std::move($1), @1.loc())); 
  }
  | keywords IDENTIFIER {
    $$ = std::move($1);
    $$.push_back(PTParsedIdentifier(std::move($2), @2.loc())); 
  }

command
  : IDENTIFIER {
    $$ = std::move(PTCommand(@1.loc(), $1));
  }
  | IDENTIFIER keywords {
    $$ = std::move(PTCommand(@1.loc(), $1));
    $$.set(std::move($2));
  }
  | IDENTIFIER LPAREN exprlist RPAREN {
    $$ = std::move(PTCommand(@1.loc(), $1));
    $$.set(std::move($3));
  }
  | IDENTIFIER keywords LPAREN exprlist RPAREN {
    $$ = std::move(PTCommand(@1.loc(), $1));
    $$.set(std::move($2));
    $$.set(std::move($4));
  }
  | IDENTIFIER opt_broken_parameter_list {
    $$ = std::move(PTCommand(@1.loc(), $1));
    $$.set(std::move($2.params));
  }
  | IDENTIFIER keywords opt_broken_parameter_list {
    $$ = std::move(PTCommand(@1.loc(), $1));
    $$.set(std::move($2));
    $$.set(std::move($3.params));
  }
  | IDENTIFIER LPAREN exprlist RPAREN opt_broken_parameter_list {
    $$ = std::move(PTCommand(@1.loc(), $1));
    $$.set(std::move($3));
    $$.set(std::move($5.params));
  }
  | IDENTIFIER keywords LPAREN exprlist RPAREN opt_broken_parameter_list {
    $$ = std::move(PTCommand(@1.loc(), $1));
    $$.set(std::move($2));
    $$.set(std::move($4));
    $$.set(std::move($6.params));
  }

control_block_build
  : CONTROL NEWLINE {
  }
  | control_block_build NEWLINE {
  }
  | control_block_build SAVE saves NEWLINE {
    // This has to be defined separately because 
    // the syntax of save command is different 
    // from the rest of commands. 
    auto cmd = PTCommand(@2.loc(), saveCmd);
    PTSaves s;
    s.add(std::move($3));
    cmd.set(std::move(s));
    tables.addCommand(std::move(cmd));
  }
  | control_block_build analysis NEWLINE {
    // Analysis also has a special syntax. 
    tables.addCommand(std::move($2));
  } 
  | control_block_build command NEWLINE {
    tables.addCommand(std::move($2));
  }
  
control_block
  : control_block_build ENDC NEWLINE {
  }

%%

namespace NAMESPACE::dflparse {

// Error reporting
void Parser::error( const Parser::location_type &l, const std::string &err_message ) {
   status.set(Status::Syntax, ("Parser "+err_message));
   status.extend(l.loc());
} 

}
