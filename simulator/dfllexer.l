%top{
/* MSVC does not define __STDC_VERSION__ which causes integer limit macros 
   to be defined by flex. These macros are later redefined by the included 
   header files and cause warnings. 
*/
#ifdef SIMMSC
#define __STDC_VERSION__ 199901L
#endif
}

%{
/* Header files */
#include <string>
#include <cmath>

/* Implementation of yyFlexScanner. Specify the prototype of lexer call. */ 
#include "dflscanner.h"
#undef  YY_DECL
#define YY_DECL int dflparse::Scanner::yylex(dflparse::Parser::semantic_type * const lval, dflparse::Parser::location_type *loc )

/* typedef to make the returns for the tokens shorter */
using token = dflparse::Parser::token;

/* define yyterminate as this instead of NULL */
#define yyterminate() return( token::END )

/* msvc2010 requires that we exclude this header file. */
#define YY_NO_UNISTD_H

/* This is called before every user action. It updates location. */
#define YY_USER_ACTION loc->step(); loc->columns(yyleng); 

/* This is called at lexer consruction- */
#define YY_USER_INIT loc->begin.setFileStack(tables.fileStack(), fileStackPosition);
%}


/* Enable start condition stack */
%option stack

/* #define FLEX_DEBUG, turn on code for debugging the lexer */
%option debug

/* Disable default rule, need to match everything manually */
%option nodefault

/* Lexer class name */
%option yyclass="Scanner"

/* No yywrap() defined. yywrap() is called at end of input. */
%option noyywrap

/* C++ lexer */
%option c++

/* Do not call isatty() */
%option never-interactive       

/* Exclusive conditions */
%x LINESTART BODY QUOTED LONGQUOTED INC INCEND LIBSECTION LIBEND SECLOOK SECLOOKNAME SECLOOKEND LONGCOMMENT

/* State transitions
-> = push
=> = pop 
.. = BEGIN

// include directive
// comments are allowed only at the end
LINESTART -> INC -> QUOTED .. INCEND => INC => LINESTART
             got 'include', must push state so that QUOTED can check the previous state
                    gout double quote
                              got include file name, looking for 'section'
                                               section not given, switched to include file
// include directive with section
LINESTART -> INC -> QUOTED .. INCEND .. LIBEND .. SECLOOK ..
                                        got section name
                                                  switched to lib file, look for section
.. SECLOOK .. SECLOOKNAME .. SECLOOKEND => INC => LINESTART (found section)
              got 'section', look for name
                             got section name, look for end of line
                                                  switched to library section
include "<filename>":
    LINESTART 
        found include: ->
    INC 
        spaces: skip
        newline: error
        found double quotes: -> 
        everything else: unmatched error
    QUOTED
        process quoted string
        found closing double quotes: ..
    INCEND
        newline: => INC =>, switch to include file
        section<optional spaces>=: .. LIBSECTION
        spaces: skip
        //: eat up comment up to newline
        /*: eat up comment up to end of comment
        everything else: unmatched error
    LINESTART

include "filename" section<optional spaces>=<optional spaces><identifier>
after INCEND we have
    LIBSECTION
        spaces: skip
        identifier: collect, ..
        newline: unexpected error
        everything else: unmatched error
    LIBEND
        newline: ..
        spaces: skip
        //: eat up comment up to newline
        /*: eat up comment up to end of comment
        everything else: unmatched error
    SECLOOK
        <begin of line><optional spaces>section<optional spaces>: ..
        newline: count lines
        everything else: skip
    SECLOOKNAME
        identifier: ..
        newline: unexpected error
        everything else: unmatched error
    SECLOOKEND
        newline: => INC =>, switch to include file
        spaces: skip
        //: eat up comment up to newline
        /*: eat up comment up to end of comment
        everything else: unmatched error
    LINESTART
*/

%%
%{          /** Code executed at the beginning of yylex **/
            // Copy matched text to parser semantic value
            yylval = lval;

            // Uncomment to debug lexer
            // yy_flex_debug = 1;
            loc->begin.setFileStack(tables.fileStack(), fileStackPosition);
            loc->end.setFileStack(tables.fileStack(), fileStackPosition);

            if (atBeginning) {
                atBeginning = false;
                if (inputType == InputNetlist) {
                    return(token::INNETLIST);
                } else if (inputType == InputExpression) {
                    return(token::INEXPR);
                }
            }
%}

.*          {
                yylval->build<std::string>(yytext);
                return(token::TITLE);
            }
\n          {
                loc->lines();
                BEGIN(LINESTART);
                return(token::NEWLINE);
            }
<LINESTART,BODY,QUOTED,LONGQUOTED,INC,INCEND,LIBSECTION,LIBEND,SECLOOK,SECLOOKNAME,SECLOOKEND,LONGCOMMENT>\r {
                // Eat up \r characters
            }
<LINESTART,BODY>\/\/[^\n]* {} // Eat up comment, do not include newline
<LINESTART,BODY>\/\* { // Eat up long comment
                yy_push_state(LONGCOMMENT);
            }
<LINESTART,BODY>[ \t]* {} // Skip spaces
<LONGCOMMENT>. {} // Eat ordinary characters in long comment
<LONGCOMMENT>\n { // Count lines in long comment
                loc->lines();
            }
<LONGCOMMENT>\*\/ { // End of long comment, go back to previous state
                yy_pop_state();
            }
<LONGQUOTED>\>\>\>[a-zA-Z0-9_]+ {
                // Possible string end
                int a=1;
                if (marker == yytext+3) {
                    yy_pop_state();
                    loc->begin = stringStart;
                    yylval->build<std::string>(std::string(sbuf));
                    return(token::STRING);
                } else {
                    // Incorrect marker
                    sbuf += yytext;
                }
            }
<LONGQUOTED>\>{1,3} {
                // 1-3 '>'
                sbuf += yytext;
            }
<LONGQUOTED>[^\>\n]+ {
                // Not '>' and not newline
                sbuf += yytext;
            }
<LONGQUOTED>\n {
                // Newline
                sbuf += yytext;
                loc->lines();
            }
<QUOTED>\"  {
                // String end (reached a double quote)
                if (yy_top_state()==INC) {
                    // If the top state on the stack is INC go to INCEND state
                    BEGIN(INCEND);
                } else {
                    yy_pop_state();
                    loc->begin = stringStart;
                    yylval->build<std::string>(std::string(sbuf));
                    return(token::STRING);
                }
            }
<QUOTED>\n  {   // Newline in quoted string
                error(*loc, "Syntax error, unexpected newline in string literal.");
                return(token::YYerror);
            }
<QUOTED>\\[0-7]{1,3} { // Octal sequence in quoted string (1-3 digits)
                int result = std::stoul(yytext+1, nullptr, 8);
                if (result > 0xff) {
                    yy_pop_state();
                    error(*loc, "Syntax error, bad octal sequence in string literal.");
                    return(token::YYerror);
                }
                sbuf += ((char)result);
            }
<QUOTED>\\[0-9]+ { // Octal sequence in quoted string (everything else)
                yy_pop_state();
                error(*loc, "Syntax error, bad escape sequence in string literal.");
                return(token::YYerror);
            }
<QUOTED>\\n { sbuf += "\n"; }
<QUOTED>\\t { sbuf += "\t"; }
<QUOTED>\\r { sbuf += "\r"; }
<QUOTED>\\b { sbuf += "\b"; }
<QUOTED>\\f { sbuf += "\f"; }
<QUOTED>[^\\\n\"]+ { // Anything, but backslash, newline, and quote in quoted string
              sbuf += yytext;
            }
<QUOTED>\\. { // Other escaped non-newline characters in quoted string
                sbuf += yytext[1];
            }
<QUOTED>\\\n { // Escaped newline in quoted string
                loc->lines();
            }
<QUOTED>.   { // Everything else in double quotes
                error(*loc, "Syntax error, unexpected string in string literal.");
                return(token::YYerror);
            }
<LINESTART>include {
                // Go to include mode
                yy_push_state(INC);
            }
<INC>[ \t]* {} // Skip spaces after include
<INC>\"     {
                // Found quote in include mode
                sbuf = "";
                stringStart = loc->begin;
                // Read quoted string (go to QUOTED state)
                yy_push_state(QUOTED);
            }
<INC>.      {
                // Unexpected input in include mode (any character)
                error(*loc, "Syntax error, unexpected string in include directive, filename must be quoted.");
                return(token::YYerror);
            }
<INC>\n     {
                // Unexpected input in include mode (newline)
                error(*loc, "Syntax error, expected quoted filename.");
                return(token::YYerror);
            }
<INCEND>[ \t]* {} // Skip spaces after quoted string
<INCEND>section[ \t]*\= {  // Found section=
                // Go to LIBSECTION state
                BEGIN(LIBSECTION);
            }
<INCEND>\/\/[^\n]* {} // Eat up comment (//) after quoted string, do not include newline
<INCEND>\/\* { // Eat up long comment after quoted string
                // Long comment
                yy_push_state(LONGCOMMENT);
            }
<INCEND>\n { // End of include directive without section name
                // INCEND is the continuation of QUOTED that was pushed from INC
                // Update line number
                loc->lines();
                // Go back to INC state
                yy_pop_state();
                // Go back to LINESTART state
                yy_pop_state();
                
                // Switch to include file
                auto pos = tables.fileStack().addFile(std::string(sbuf), Simulator::includePath(), fileStackPosition, loc->begin.line);
                if (pos!=FileStack::badFileId) {
                    fileStackPosition = pos;
                    
                    if (Simulator::fileDebug()) {
                        Simulator::dbg() << "Opening include file '" << tables.fileStack().canonicalName(fileStackPosition) << "'.\n";
                    }
                    
                    auto& fname = tables.fileStack().canonicalName(fileStackPosition);
                    auto streamPtr = pushStream(fname, *loc);
                    if (!streamPtr) {
                        error(*loc, std::string("Failed to open include file '")+fname+"'.");
                        return(token::YYerror);
                    }
                    yypush_buffer_state(
                        yy_create_buffer(
                            *streamPtr, YY_BUF_SIZE
                        )
                    );
                    loc->end.initialize(nullptr, 1, 1, 0);
                    loc->end.setFileStack(tables.fileStack(), pos);
                    loc->begin = loc->end;
                } else {
                    loc->begin = stringStart;
                    error(*loc, "File not found.");
                    return(token::YYerror);
                }
            }
<INCEND>.   {
                // Unexpected input in include end mode 
                error(*loc, "Syntax error, unexpected string in include directive, filename must be quoted.");
                return(token::YYerror);
            }
<LIBSECTION>[ \t]* {} // Skip spaces after section=
<LIBSECTION>[a-zA-Z_$][0-9a-zA-Z_$]* { // Found an identifier
                // Save section name
                section = yytext;
                // Go to LIBEND state
                BEGIN(LIBEND);
            }
<LIBSECTION>\n {
                // Unexpected input in section name lookup
                error(*loc, "Syntax error, expected section name.");
                return(token::YYerror);
            }
<LIBSECTION>. {
                // Unexpected input in section name lookup
                error(*loc, "Syntax error, unexpected string in include directive.");
                return(token::YYerror);
            }
<LIBEND>[ \t]* {} // Skip spaces after include section name
<LIBEND>\/\/[^\n]* {} // Eat up comment (//), do not include newline
<LIBEND>\/\* { // Eat up long comment after quoted string
                // Long comment
                yy_push_state(LONGCOMMENT);
            }
<LIBEND>\n  { // End of include directive with section name
                // Switch to lib file
                // Update line number
                loc->lines();
                // Switch to library file
                auto pos = tables.fileStack().addFile(std::string(sbuf), Simulator::includePath(), fileStackPosition, loc->begin.line);
                if (pos!=FileStack::badFileId) {
                    fileStackPosition = pos;
                    tables.fileStack().setSection(pos, section);

                    if (Simulator::fileDebug()) {
                        Simulator::dbg() << "Opening library file '" << tables.fileStack().canonicalName(fileStackPosition) 
                            << "', section '" << section << "'.\n";
                    }
                    
                    auto& fname = tables.fileStack().canonicalName(fileStackPosition);
                    auto streamPtr = pushStream(fname, *loc);
                    if (!streamPtr) {
                        error(*loc, std::string("Failed to open library file '")+fname+"'.");
                        return(token::YYerror);
                    }
                    yypush_buffer_state(yy_create_buffer(
                            *streamPtr, 
                            YY_BUF_SIZE
                        )
                    );
                    loc->end.initialize(nullptr, 1, 1, 0);
                    loc->end.setFileStack(tables.fileStack(), pos);
                    loc->begin = loc->end;
                    // LIBEND is a continuation of LIB state
                    // Go to SECLOOK state
                    BEGIN(SECLOOK);
                } else {
                    loc->begin = stringStart;
                    error(*loc, "File not found.");
                    return(token::YYerror);
                }
            }
<LIBEND>.   { 
                // Unexpected input in include end mode 
                error(*loc, "Syntax error, unexpected string in library directive.");
                return(token::YYerror);
            }
<SECLOOK>^[ \t]*\section[ \t]* {
                // Found section, skipped trailing spaces
                // Go to SECLOOKNAME
                BEGIN(SECLOOKNAME);
            }
<SECLOOK>.  {}  // Skip other characters (not newline) when searching for section
<SECLOOK>\n {
                // Newline is also skipped, line counter is increased
                loc->lines();
            }
<SECLOOKNAME>[a-zA-Z_$][0-9a-zA-Z_$]* {
                // Found identifier (section name)
                if (section==yytext) {
                    // Matches section we are looking for
                    // Go to SECLOOKEND state
                    BEGIN(SECLOOKEND);
                } else {
                    // No match, go to SECLOOK state
                    BEGIN(SECLOOK);
                }
            }
<SECLOOKNAME>\n {
                // Section search, unexpected input when looking for section name
                error(*loc, "Syntax error, expected section name.");
                return(token::YYerror);
            }
<SECLOOKNAME>. {
                // Section search, unexpected input when looking for section name
                error(*loc, "Syntax error, unmatched string in section directive.");
                return(token::YYerror);
            }
<SECLOOKEND>[ \t]* { // Section found, skip whitespace after section name until and including newline
            } 
<SECLOOKEND>[ \t]*\/\/[^\n]* { // Section found, eat up comment, do not include newline
            } 
<SECLOOKEND>[ \t]*\/\* { // Section found, eat up long comment
                // Long comment
                yy_push_state(LONGCOMMENT);
            }
<SECLOOKEND>\n {
                // Update line number
                loc->lines();
                // SECLOOKEND is a continuation of QUOTED state that was pushed from LIB state
                // Go back to LIB state
                yy_pop_state();
                // Go back to LINESTART state
                yy_pop_state();
            }
<SECLOOKEND>. {
                // Unexpected input in section found mode 
                error(*loc, "Syntax error, unmatched string in section directive.");
                return(token::YYerror);
            }
<LINESTART>endsection {
                // Are we in a file section
                if (tables.fileStack().sectionName(fileStackPosition).size()>0) {
                    // Yes, finished with file
                    // We are in LINESTART state, which is where we must also be when we 
                    // return to calling file because we departed after a newline. 
                    // Restore parent buffer
                    yypop_buffer_state();
                    // Restore parent location, close stream
                    *loc = popStream();
                    // Go back in file stack
                    fileStackPosition = tables.fileStack().parentId(fileStackPosition);
                } else {
                    error(*loc, "Syntax error, endsection directive reached in a file not included with a library directive.");
                    return(token::YYerror);
                }
            }
<LINESTART>section[ \t]* {
                // Found section, skipped trailing spaces
                error(*loc, "Syntax error, unexpected section directive. You probably incldued a library file.");
                return(token::YYerror);
            }

<BODY>\(    {
                inParen++;
                return(token::LPAREN);
            }

<BODY>\)    {
                inParen--;
                return(token::RPAREN);
            }

<BODY>\[    {
                inBracket++;
                return(token::LBRACKET);
            }
<BODY>\]    {
                inBracket--;
                return(token::RBRACKET);
            }

<BODY>\=    {
                return(token::ASSIGN);
            }
<BODY>\,    {
                return(token::COMMA);
            }

<BODY>0|([1-9][0-9]*) {
                yylval->build<Int>(std::stoi(yytext, nullptr, 10));
                return(token::INTEGER);
            }

<BODY>0[xX][0-9]* {
                yylval->build<Int>(std::stoul(yytext, nullptr, 16));
                return(token::HEXINTEGER);
            }

<BODY>(0|([1-9][0-9]*))?\.([0-9]*)? {
                // xxx.xxx
                // xxx.
                yylval->build<double>(std::stod(yytext));
                return(token::FLOAT);
            }
<BODY>(0|([1-9][0-9]*))?(\.([0-9]*)?)?([eE][+-]?[0-9]+) {
                // [.xxx | xxx[.[xxx]] | xxx] e [+|-] xxx
                yylval->build<double>(std::stod(yytext));
                return(token::FLOAT);
            }
<BODY>(\.[0-9]*)|((0|([1-9][0-9]*))(\.([0-9]*)?)?)([munpfakKMGT][a-zA-Z_]*) { 
                // .xxx
                // xxx[.[xxx]]
                // xxx
                // with SI-prefix and optional unit
                double s;
                int i;
                for(i=0; yytext[i]; i++) {
                    if (!(yytext[i]>='0' && yytext[i]<='9') && yytext[i]!='.') {
                        break;
                    }
                }
                // i is the index of SI unit prefix
                switch (yytext[i]) {
                    case 'm': s=1e-3; break;
                    case 'u': s=1e-6; break;
                    case 'n': s=1e-9; break;
                    case 'p': s=1e-12; break;
                    case 'f': s=1e-15; break;
                    case 'a': s=1e-18; break;
                    case 'k': s=1e3; break;
                    case 'K': s=1e3; break;
                    case 'M': s=1e6; break;
                    case 'G': s=1e9; break;
                    case 'T': s=1e12; break;
                    default: 
                        error(*loc, "Syntax error, incorrectly formed floating point number.");
                        return(token::YYerror);
                }
                yylval->build<double>(std::stod(std::string(yytext, 0, i))*s);
                return(token::FLOAT);
            }
<BODY>(?i:inf)|(?i:nan) {
                yylval->build<double>(std::stod(yytext));
                return(token::FLOAT);
            }
<BODY>\<\<\<[a-zA-Z0-9_]+[\r\t ]*\n {
                loc->lines();
                sbuf = "";
                marker = yytext+3;
                // Trim trailing spaces
                auto pos = marker.find_last_not_of("\n\r\t ");
                marker = marker.substr(0, pos+1);
                stringStart = loc->begin;
                yy_push_state(LONGQUOTED);
            }
<BODY>\" {
                sbuf = "";
                stringStart = loc->begin;
                yy_push_state(QUOTED);
            }
<BODY>\+    {
                return (token::PLUS);
            }
<BODY>\-    {
                return (token::MINUS);
            }
<BODY>\*    {
                return (token::TIMES);
            }
<BODY>\/    {
                return (token::DIVIDE);
            }
<BODY>\*\*  {
                return (token::POWER);
            }
<BODY>==  {
                return (token::EQUAL);
            }
<BODY>!=  {
                return (token::NOTEQUAL);
            }
<BODY>\<=  {
                return (token::LESSEQ);
            }
<BODY>\<  {
                return (token::LESS);
            }
<BODY>\>=  {
                return (token::GREATEREQ);
            }
<BODY>\>  {
                return (token::GREATER);
            }
<BODY>\&  {
                return (token::BITAND);
            }
<BODY>\|  {
                return (token::BITOR);
            }
<BODY>\^  {
                return (token::BITEXOR);
            }
<BODY>\~  {
                return (token::BITNOT);
            }
<BODY>\>\>  {
                return (token::BITSHIFTR);
            }
<BODY>\<\<  {
                return (token::BITSHIFTL);
            }
<BODY>\&\&  {
                return (token::AND);
            }
<BODY>\|\|  {
                return (token::OR);
            }
<BODY>\!    {
                return (token::NOT);
            }
<BODY>-\>   {
                return (token::RIGHTARROW);
            }
<BODY>:     { 
                return (token::COLON);
            }
<BODY>\?    { 
                return (token::QUESTION);
            }
<BODY>;     {
                return (token::SEMICOLON);
            }

<LINESTART>model {
                BEGIN(BODY);
                return (token::MODEL);
            }
<LINESTART>global {
                // Do not switch to BODY yet. 
                BEGIN(BODY);
                return (token::GLOBAL);
            }
<LINESTART>ground {
                BEGIN(BODY);
                return (token::GROUND);
            }
<LINESTART>analysis {
                BEGIN(BODY);
                return (token::ANALYSIS);
            }
<LINESTART>sweep {
                BEGIN(BODY);
                return (token::SWEEP);
            }
<LINESTART>load {
                BEGIN(BODY);
                return (token::LOAD);
            }
<LINESTART>embed {
                BEGIN(BODY);
                return (token::EMBED);
            }
<LINESTART>subckt {
                BEGIN(BODY);
                return (token::SUBCKT);
            }

<LINESTART>ends {
                BEGIN(BODY);
                return (token::ENDS);
            }
<LINESTART>parameters {
                BEGIN(BODY);
                return (token::PARAMETERS);
            }
<LINESTART>control {
                BEGIN(BODY);
                inControl = true;
                return (token::CONTROL);
            }
<LINESTART>endc {
                BEGIN(BODY);
                inControl = false;
                return (token::ENDC);
            }
<LINESTART>save {
                BEGIN(BODY);
                if (inControl) {
                    return (token::SAVE);
                } else {
                    yylval->build<Id>(yytext);
                    return (token::IDENTIFIER);
                }
            }
<LINESTART,BODY>[a-zA-Z_$][0-9a-zA-Z_$]* { 
                // Identifier
                BEGIN(BODY);
                yylval->build<Id>(yytext);
                return(token::IDENTIFIER);
            }
<LINESTART,BODY>'[0-9a-zA-Z_$:\.\(\)]+' { 
                // Identifier in single quotes, any nonempty alphanumeric sequence with _, $, :, ., (, )
                BEGIN(BODY);
                auto n = std::strlen(yytext);
                yylval->build<Id>(std::string(yytext+1, n-2));
                return(token::IDENTIFIER);
            }

<LINESTART,BODY>\n { // Newline, return token
                // Update line number
                loc->lines();
                if (inParen<=0 && inBracket<=0) {
                    BEGIN(LINESTART);
                    return( token::NEWLINE );
                }
            }
<LINESTART,BODY>\\\n { // Escaped newline, no token
                // Update line number, no token returned because of line continuation
                loc->lines();
            }

<LINESTART,BODY>. { 
                error(*loc, "Syntax error, unexpected string.");
                return(token::YYerror);
            }

<LINESTART,BODY><<EOF>> {
                if (fileStackPosition>0) {
                    const std::string& tmp = tables.fileStack().sectionName(fileStackPosition);
                    if (tmp.size()>0) {
                        // Error, EOF within a section
                        loc->begin = loc->end;
                        error(*loc, std::string("Syntax error, unexpected end of file in section '")+tmp+"'.");
                        return(token::YYerror);
                    } else {
                        // EOF in included file
                        // Restore parent buffer
                        yypop_buffer_state();
                        // Restore parent location, close stream
                        *loc = popStream();
                        // Go back in file stack
                        fileStackPosition = tables.fileStack().parentId(fileStackPosition);
                    }
                } else {
                    // EOF in top level file
                    loc->begin = loc->end;
                    return(token::END);
                } 
            }            
<INC,INCEND,LIBSECTION,LIBEND,LONGQUOTED,QUOTED,LONGCOMMENT><<EOF>> {
                loc->begin = loc->end;
                error(*loc, "Syntax error, unexpected end of file.");
                return(token::YYerror);
            }
<SECLOOK,SECLOOKNAME,SECLOOKEND><<EOF>> {
                const std::string& tmp = tables.fileStack().sectionName(fileStackPosition);
                loc->begin = loc->end;
                error(*loc, std::string("Syntax error, unexpected end of file while looking for section '")+tmp+"'.");
                return(token::YYerror);
            }
%%
