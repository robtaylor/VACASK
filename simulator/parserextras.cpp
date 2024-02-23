#include "parserextras.h"
#include "common.h"


namespace NAMESPACE {

PTEmbed::PTEmbed() {
}

PTEmbed::PTEmbed(const Loc& l, std::string&& filename, std::string&& contents) 
    : loc_(l), filename_(std::move(filename)), contents_(std::move(contents)) {
}

PTEmbed::~PTEmbed() {
}

std::ostream& operator<<(std::ostream& os, const PTEmbed& e) {
    os << "embed \"" << e.filename_ << "\" <<<FILE\n" << e.contents_ << ">>>FILE";
    return os;
}


PTCommand::PTCommand() {
}

PTCommand::PTCommand(const Loc& l, Id name) 
    : loc_(l), name_(name) {
}

void PTCommand::add(PTIdentifierList&& kw) {
    keywords_ = std::move(kw);
}

void PTCommand::add(PTParameters&& args) {
    args_ = std::move(args);
}

void PTCommand::add(std::vector<Rpn>&& exprs) {
    expressions_ = std::move(exprs);
}

void PTCommand::add(PTSaves&& s) {
    saves_ = std::move(s);
}

PTCommand::~PTCommand() {
}

std::ostream& operator<<(std::ostream& os, const PTCommand& c) {
    os << c.name_;
    if (c.keywords_.size()>0) {
        os << " " << c.keywords_;
    }
    if (c.saves_.saves().size()==0) {
        if (c.expressions_.size()>0) {
            os << " (";
            bool first = true;
            for(auto& it : c.expressions_) {
                if (!first) {
                    os << ", ";
                }
                os << it;
                first = false;
            }
            os << ")";
        }
        os << " " << c.args_;
    } else {
        os << " " << c.saves_;
    }
    return os;
}


ParserExtras::ParserExtras() {
}

ParserExtras::~ParserExtras() {
}

void ParserExtras::addCommand(PTAnalysis&& a) {
    control_.push_back(std::move(a));
}

void ParserExtras::addCommand(PTCommand&& c) {
    control_.push_back(std::move(c));
}

void ParserExtras::addEmbed(PTEmbed&& e) {
    embed_.push_back(std::move(e));
}

bool ParserExtras::verify(Status& s) const {
    // Check duplicate analyses
    std::unordered_map<Id,const PTAnalysis*> anmap;
    for(auto& it : control_) {
        if (!std::holds_alternative<PTAnalysis>(it)) {
            continue;
        }
        auto& an = std::get<PTAnalysis>(it);
        auto [itEx, inserted] = anmap.insert({an.name(), &an});
        if (!inserted) {
            s.set(Status::Redefinition, "Analysis '"+std::string(an.name())+"' redefinition.");
            s.extend(an.location());
            s.extend("Analysis was first defined here");
            s.extend(itEx->second->location());
            return false;
        }
    }

    return true;
}

void ParserExtras::dump(int indent, std::ostream& os) {
    std::string pfx = std::string(indent, ' ');
    if (embed_.size()>0) {
        for(auto& it : embed_) {
            os << pfx << it << "\n";
        }
        os << "\n";
    }

    if (control_.size()>0) {
        os << pfx << "control\n";
        for(auto& it : control_) {
            if (std::holds_alternative<PTAnalysis>(it)) {
                std::get<PTAnalysis>(it).dump(indent+2, os); 
            } else {
                os << pfx << "  " << std::get<PTCommand>(it) << "\n";
            }
        }
        os << pfx << "endc\n";
    }
}

}
