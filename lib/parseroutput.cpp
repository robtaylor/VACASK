#include "parseroutput.h"
#include "simulator.h"
#include "status.h"
#include "common.h"


namespace NAMESPACE {

std::ostream& operator<<(std::ostream& os, const PTParsedIdentifier& obj) {
    os << std::string(obj.id_);
    return os;
}


std::ostream& operator<<(std::ostream& os, const PTIdentifierList& obj) {
    for(auto it=obj.cbegin(); it!=obj.cend(); ++it) {
        if (it!=obj.cbegin()) 
            os << " ";
        if (it->name()!=Id::none) {
            os << it->name();
        } else {
            os << "<null>";
        }
    }
    return os;
}


void PTParameterValue::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << std::string(id_) << "=" << val_;
}


void PTParameterExpression::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << std::string(id_) << "=" << rpn_.str();
}


std::ostream& operator<<(std::ostream& os, const PTParameters& obj) {
    for(auto it=obj.values_.cbegin(); it!=obj.values_.cend(); ++it) {
        os << (it->name()) << "=" << it->val() << " " ;
    }
    for(auto it=obj.expressions_.begin(); it!=obj.expressions_.end(); ++it) {
        os << (it->name()) << "=" << it->rpn() << " " ;
    }
    
    return os;
}


void PTModel::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "model " << (modelName_) << " " << (deviceName_) << " ";
    os << parameters_ << "\n";
}


void PTInstance::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << (instanceName_) << " (" << connections_ << ") ";
    os << masterName_ << " " << parameters_ << "\n";
}


void PTBlockSequence::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    bool first = true;
    for(auto& entry : entries_) {
        auto& [loc, rpn, block] = entry;
        if (rpn.size()>0) {
            if (first) {
                os << pfx << "@if " << rpn << "\n";
            } else {
                os << pfx << "@elseif " << rpn << "\n";
            }
        } else {
            os << pfx << "@else\n";
        }
        block.dump(indent+2, os);
        first = false;
    }
    os << pfx << "@end\n";
}


void PTBlock::dump(int indent, std::ostream& os) const {
    for(auto& mod : models_) {
        mod.dump(indent, os);
    }
    for(auto& inst : instances_) {
        inst.dump(indent, os);
    }
    if (hasBlockSequences()) {
        for(auto& seq: *blockSequences_) {
            seq.dump(indent, os);
        }
    }
}

void PTSubcircuitDefinition::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');

    bool isToplevel = !modelName_;

    if (!isToplevel) {
        os << pfx << "subckt " << (modelName_) << " (" << terminals_ << ")\n";
    }
    
    if (parameters_.count()>0) {
        os << pfx << (isToplevel ? "" : "  ") << "parameters " << parameters_ << "\n";
    }
    if (subDefs_.size()>0) {
        os << "\n";
    }
    for(auto it=subDefs_.begin(); it!=subDefs_.end(); ++it) {
        it->get()->dump(isToplevel ? indent : indent+2, os);
    }
    if (root_.models().size()>0) {
        os << "\n";
    }
    root_.dump(isToplevel ? indent : indent+2, os);
    
    if (!isToplevel) {
        os << pfx << "ends\n";
    }
    
    os << "\n";
}


void PTLoad::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "load \"" << file_ << "\"";
    os << " " << (parameters_);
    os << "\n";
}


std::ostream& operator<<(std::ostream& os, const PTSave& s) {
    os << std::string(s.typeName()) << "(";
    if (s.id[0]) {
        os << "\"" << std::string(s.id[0]) << "\"";
    }
    if (s.id[1]) {
        os << "," << "\"" << std::string(s.id[1]) << "\"";
    }
    os << ")";
    
    return os;
}


std::ostream& operator<<(std::ostream& os, const PTSaves& s) {
    if (s.saves_.size()>0) {
        for(auto it=s.saves_.cbegin(); it!=s.saves_.cend(); ++it) {
            os << *it << " ";
        }
    }
    
    return os;
}


std::ostream& operator<<(std::ostream& os, const PTSweep& s) {
    os << "sweep " << std::string(s.name_) << " " << s.parameters_;
    
    return os;
}


void PTSweeps::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    for(auto it=sweeps().cbegin(); it!=sweeps().cend(); ++it) {
        os << pfx << *it << "\n";
    }
}


void PTAnalysis::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    if (sweeps_.size()>0) {
        for(auto it=sweeps_.cbegin(); it!=sweeps_.cend(); ++it) {
            os << pfx << *it << "\n";
        }
    }
    os << pfx << (sweeps_.size()>0 ? "  " : "");
    os << "analysis " << std::string(name_) << " " << std::string(typeName_) << " ";
    os << parameters_ << "\n";
}


std::ostream& operator<<(std::ostream& os, const PTEmbed& e) {
    os << "embed \"" << e.filename_ << "\" <<<FILE\n" << e.contents_ << ">>>FILE";
    return os;
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


bool PTParameters::verify(int level, Status& s) const {
    std::unordered_map<Id,Loc> puniq;
    for(auto& it : values_) {
        auto [itPrev, inserted] = puniq.insert({it.name(), it.location()});
        if (!inserted) {
            s.set(Status::Redefinition, "Parameter '"+std::string(it.name())+"' redefinition.");
            s.extend(it.location());
            if (itPrev->second) {
                s.extend("Parameter was first defined here");
                s.extend(itPrev->second);
            }
            return false;
        }
    }
    for(auto& it : expressions_) {
        auto [itPrev, inserted] = puniq.insert({it.name(), it.location()});
        if (!inserted) {
            s.set(Status::Redefinition, "Parameter '"+std::string(it.name())+"' redefinition.");
            s.extend(it.location());
            if (itPrev->second) {
                s.extend("Parameter was first defined here");
                s.extend(itPrev->second);
            }
            return false;
        }
    }
    return true;
}

bool PTModel::verify(int level, Status& s) const {
    if (!parameters_.verify(level, s)) {
        if (!loc) {
            s.extend("  in model '"+std::string(modelName_)+"'");
        }
        return false;
    }
    return true;
}

bool PTInstance::verify(int level, Status& s) const {
    if (!parameters_.verify(level, s)) {
        if (!loc) {
            s.extend("  in instance '"+std::string(instanceName_)+"'");
        }
        return false;
    }
    return true;
}

bool PTBlock::verify(int level, Status& s) const {
    // Check models
    for(auto& mod : models_) {
        if (!mod.verify(level, s)) {
            return false;
        }
    }

    // Check instances
    for(auto& inst : instances_) {
        if (!inst.verify(level, s)) {
            return false;
        }
    }

    // Recurse into all blocks in all block sequences
    if (blockSequences_) {
        int cnt=1;
        for(auto& blkSeq : *blockSequences_) {
            if (!blkSeq.verify(level, s)) {
                if (!loc) {
                    s.extend("  in block sequence '"+std::to_string(cnt)+"'");
                }
                return false;
            }
            cnt++;
        }
    }
    return true;
}

bool PTBlockSequence::verify(int level, Status& s) {
    // Check all blocks in all block sequences
    int cnt = 1;
    for(auto& blkSeqEntry : entries_) {
        auto& [loc, cond, blk] = blkSeqEntry;
        if (!blk.verify(level, s)) {
            if (!loc) {
                s.extend("  in block '"+std::to_string(cnt)+"'");
            }
            return false;
        }
        cnt++;
    }
    return true;
}

bool PTSubcircuitDefinition::verify(int level, Status& s) const {
    // Subcircuit terminal can have the same name as a global node. 
    // In that case it is simply a local node name. It does not 
    // behave as a global node. 
    // So do not check against global nodes. 

    // Check for duplicate terminal names
    std::unordered_set<Id> tset;
    for(auto& term : terminals_) {
        auto [exIt, inserted] = tset.insert(term.name());
        if (!inserted) {
            s.set(Status::Conflicting, "Terminal '"+std::to_string(term.name())+"' is not unique.");
            if (loc) {
                s.extend(loc);
            } else {
                s.extend("  in subcircuit definition '"+std::string(name())+"'");
            }
            
            return false;
        }
    }

    // Do not check name conflicts between instances, models, and subcircuits
    // Elaboration will handle that

    // Check parameters
    if (!parameters_.verify(level, s)) {
        if (!loc) {
            s.extend("  in subcircuit definition '"+std::string(name())+"'");
        }
        return false;
    }
    
    // Check root block
    if (!root_.verify(level, s)) {
        if (!loc) {
            s.extend("  in subcircuit definition '"+std::string(name())+"'");
        }
        return false;
    }

    // Check all subcircuit definitions within this definition
    for(auto& subDef : subDefs_) {
        if (!subDef->verify(level, s)) {
            if (!loc) {
                s.extend("  in subcircuit definition '"+std::string(name())+"'");
            }
            return false;
        }
    }
    return true; 
}

bool PTLoad::verify(int level, Status& s) const {
    // Check for parameters whose values are defined with expressions
    if (parameters_.expressionCount()>0) {
        s.set(Status::BadArguments, "Load directive parameters must be constants.");
        if (loc) {
            s.extend(loc);
        } else {
            s.extend("  in load directive for '"+file_+"'.");
        }
        return false;
    }

    // Check for parameter duplicates
    if (!parameters_.verify(level, s)) {
        if (!loc) {
            s.extend("  in load directive for '"+file_+"'.");
        }
        return false;
    }

    return true;
}

bool PTSweep::verify(int level, Status& s) const {
    // Verify parameters
    if (!parameters_.verify(level, s)) {
        if (!loc) {
            s.extend("  in sweep '"+std::string(name_)+"'.");
        }
        return false;
    }
    return true;
}

bool PTAnalysis::verify(int level, Status& s) const {
    // Verify sweeps
    for(auto& sw : sweeps_) {
        if (!sw.verify(level, s)) {
            if (!loc) {
                s.extend("  in analysis '"+std::string(name_)+"'.");
            }
            return false;
        }
    }

    // Verify parameters
    if (!parameters_.verify(level, s)) {
        if (!loc) {
            s.extend("  in analyis '"+std::string(name_)+"'.");
        }
        return false;
    }

    return true;
}

bool ParserTables::verifyWorker(int level, Status& s) const {
    // Check for duplicate analyses (all levels)
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
            if (itEx->second->location()) {
                s.extend("Analysis was first defined here");
                s.extend(itEx->second->location());
            }
            return false;
        }
    }

    // Verify load directives
    for(auto& load : loads_) {
        if (!load.verify(level, s)) {
            return false;
        }
    }

    // Verify default subcircuit definition
    if (level>0) {
        if (!defaultSubDef_.verify(level, s)) {
            return false;
        }
    }

    return true;
}


bool ParserTables::writeEmbedded(int debug, Status& s) {
    for(auto& e : embed_) {
        // Do we have a valid location of the embed directive
        bool dump = false;
        std::string originatorFile;
        if (!e.location()) {
            // Location not available, always dump
            dump = true;
        } else {
            // Get canonical path of file with the embed directive
            auto [fs, pos, line, offset] = e.location().data();
            auto directiveLocationCanonicalPath = fs->canonicalName(pos);

            // Build name of origin file - add .origin to dumped file name
            auto originFilePath = e.filename() + ".origin";
            
            // Does origin file exist, does dumped file exist
            if (
                !std::filesystem::exists(originFilePath) ||
                !std::filesystem::exists(e.filename())
            ) {
                // No, dump
                dump = true;
            } else {
                // Read origin file, get name of the file with the embed directive that
                // produced the existing dumped file
                std::ifstream file(originFilePath);
                if (!file) {
                    // Failed to read, dump
                    dump = true;
                } else {
                    // Read origin file to get actual origin
                    originatorFile = std::string(
                        std::istreambuf_iterator<char>(file),
                        std::istreambuf_iterator<char>()
                    );
                    file.close();

                    // Does the actual origin match the embed directive location
                    if (originatorFile!=directiveLocationCanonicalPath) {
                        // No, dump
                        dump = true;
                    } else {
                        // Compare last modification of originator and dumped file
                        auto refModificationTime = std::filesystem::last_write_time(originatorFile);
                        auto fileModificationTime = std::filesystem::last_write_time(e.filename());
                        if (refModificationTime>fileModificationTime) {
                            // Originator is newer, dump
                            dump = true;
                        }
                    }
                }
            }

            // No need to dump
            if (!dump) {
                continue;
            }
            
            // Create origin file
            std::ofstream originFile(originFilePath, std::ios::out);
            if (!originFile) {
                // Failure to write origin file is silently inored
            } else {
                // Dump and close
                originFile << directiveLocationCanonicalPath; 
                originFile.close();
            }

            std::ofstream file(e.filename(), std::ios::out);
            if (!file) {
                // Failure to dump is an error
                s.set(Status::CreationFailed, "Failed to write file '"+e.filename()+"'.");
                s.extend(e.location());
                return false;
            } else {
                // Dump
                file << e.contents();
                file.close();
            }
        }
    }
    return true;
}

void ParserTables::dump(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');

    os << pfx << title_ << "\n\n";

    if (loads_.size()>0) {
        for(auto it=loads_.begin(); it!=loads_.end(); ++it) {
            it->dump(indent, os);
        }
        os << "\n";
    }
    
    if (groundNodes_.size()>0) {
        os << pfx << "ground";
        for(auto it=groundNodes_.begin(); it!=groundNodes_.end(); ++it) {
            os << " " << *it;
        }
        os << "\n\n";
    }
    
    if (globalNodes_.size()>0) {
        os << pfx << "global";
        for(auto it=globalNodes_.begin(); it!=globalNodes_.end(); ++it) {
            os << " " << *it;
        }
        os << "\n\n";
    }

    defaultSubDef_.dump(indent, os);
    os << "\n";

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
