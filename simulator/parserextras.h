#ifndef __PARSEREXTRAS_DEFINED
#define __PARSEREXTRAS_DEFINED

#include <string>
#include <cstddef>
#include <istream>
#include <unordered_map>
#include "location.h"
#include "parseroutput.h"
#include "common.h"


namespace NAMESPACE {

class PTEmbed {
public:
    PTEmbed();
    PTEmbed(const Loc& l, std::string&& filename, std::string&& contents);
    ~PTEmbed();

    PTEmbed           (const PTEmbed&)  = delete;
    PTEmbed           (      PTEmbed&&) = default;
    PTEmbed& operator=(const PTEmbed&)  = delete;
    PTEmbed& operator=(      PTEmbed&&) = default;

    Loc location() const { return loc_; };
    const std::string& filename() const { return filename_; }; 
    const std::string& contents() const { return contents_; }; 

    friend std::ostream& operator<<(std::ostream& os, const PTEmbed& e);

private:
    Loc loc_;
    std::string filename_;
    std::string contents_;
};


class PTCommand {
public: 
    PTCommand();
    PTCommand(const Loc& l, Id name);
    ~PTCommand();

    PTCommand           (const PTCommand&)  = delete;
    PTCommand           (      PTCommand&&) = default;
    PTCommand& operator=(const PTCommand&)  = delete;
    PTCommand& operator=(      PTCommand&&) = default;

    void add(PTIdentifierList&& kw);
    void add(PTParameters&& args);
    void add(std::vector<Rpn>&& exprs);
    void add(PTSaves&& s);

    Loc location() const { return loc_; };
    Id name() const { return name_; };
    const PTIdentifierList& keywords() const { return keywords_; }; 
    const std::vector<Rpn>& expressions() const { return expressions_; }; 
    const PTParameters& args() const { return args_; }; 
    const PTSaves& saves() const { return saves_; };
    PTSaves& saves() { return saves_; };
    
    friend std::ostream& operator<<(std::ostream& os, const PTCommand& c);

private:
    Loc loc_;
    Id name_;
    std::vector<Rpn> expressions_;
    PTIdentifierList keywords_;
    PTParameters args_;
    PTSaves saves_;
};

typedef std::variant<PTAnalysis, PTCommand> ControlEntry;
typedef std::vector<ControlEntry> PTControl;


class ParserExtras {
public:
    ParserExtras();
    ~ParserExtras();

    ParserExtras           (const ParserExtras&)  = delete;
    ParserExtras           (      ParserExtras&&) = default;
    ParserExtras& operator=(const ParserExtras&)  = delete;
    ParserExtras& operator=(      ParserExtras&&) = default;

    void addEmbed(PTEmbed&& e);
    void addCommand(PTSaves&& s);
    void addCommand(PTCommand&& c);
    void addCommand(PTAnalysis&& a);
    
    void setExpr(Rpn&& e) { expr_ = std::move(e); };

    const PTControl& control() const { return control_; }; 
    PTControl& control() { return control_; }; 
    const std::vector<PTEmbed>& embed() const { return embed_; };
    Rpn& expr() { return expr_; };

    // Post-parse checks
    bool verify(Status& s=Status::ignore) const;

    void dump(int indent, std::ostream& os);
    
private:
    std::vector<PTEmbed> embed_;
    PTControl control_;
    Rpn expr_;
};

}

#endif
