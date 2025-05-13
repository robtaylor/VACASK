#include <cmath>
#include <exception>
#include "rpneval.h"
#include "rpnbuiltin.h"
#include "rpnfunctor.h"
#include "common.h"


namespace NAMESPACE {

bool RpnEvaluator::isConstant(const Rpn& expr) const {
    for(auto e=expr.cbegin(); e!=expr.cend(); ++e) {
        if (e->type()==Rpn::TIdentifier) {
            // Is it a constant, OK, otherwise not
            if (!ContextStack::isConstant(e->get<Rpn::Identifier>().name)) {
                return false;
            }
        } else if (e->type()==Rpn::TFunctionCall) {
            // Is it a pure function, OK, otherwise not
            auto& fcall = e->get<Rpn::FunctionCall>();
            if (!ContextStack::getBuiltin(fcall.name)->pure) {
                return false;
            }
        }
    }
    return true;
}

void RpnEvaluator::appendLocation(Status& s, const Loc& p) {
    s.extend(p);
}

std::tuple<bool, bool> RpnEvaluator::checkCondition() {
    Value* cond = stack_.get(); 
    switch (cond->type()) {
        case Value::Type::Int: return std::make_tuple(true, makeBool(cond->val<Int>()));
        case Value::Type::Real: return std::make_tuple(true, makeBool(cond->val<Real>()));
        case Value::Type::String: return std::make_tuple(true, cond->val<String>().size()>0);
        default: return std::make_tuple(false, false);
    }
}

bool RpnEvaluator::evaluate(const Rpn& rpn, Value& result, Status& s) {
    stack_.clear();
    Value *v1p, *v2p;
    // Increment manually at the end of the loop
    for(auto e=rpn.cbegin(); e!=rpn.cend(); ) {
        const Loc& loc = rpn.location(*e);
        bool isJump = false; 
        switch (e->type()) {
            case Rpn::TValue: {
                stack_.push(e->get<Value>());
                break;
            }
            case Rpn::TIdentifier: {
                auto& id = e->get<Rpn::Identifier>();
                auto* valPtr = contextStack_.get(id.name);
                if (!valPtr) {
                    s.set(Status::NotFound, std::string("Variable or constant '")+std::string(id.name)+"' not defined.");
                    s.extend(loc);
                    return false;
                }
                stack_.push(*valPtr);
                break;
            }
            case Rpn::TOp: try {
                    // Operator
                    auto& op = e->get<Rpn::Op>();
                    bool coreSt = false;
                    switch (op.code) {
                        case Rpn::OpUMinus:
                            coreSt = mathFuncComp1<FwUminus>(stack(), 1, s); break;
                        case Rpn::OpPlus:
                            coreSt = mathFuncComp2<FwPlus>(stack(), 2, s); break;
                        case Rpn::OpMinus:
                            coreSt = mathFuncComp2<FwMinus>(stack(), 2, s); break;
                        case Rpn::OpTimes:
                            coreSt = mathFuncComp2<FwTimes>(stack(), 2, s); break;
                        case Rpn::OpDivide:
                            coreSt = mathFuncComp2<FwDivide>(stack(), 2, s); break;
                        case Rpn::OpPower:
                            coreSt = mathFuncComp2<FwPower>(stack(), 2, s); break;
                        case Rpn::OpEqual:
                            coreSt = mathRelOpComp2<FwEqual>(stack(), 2, s); break;
                        case Rpn::OpNotEqual:
                            coreSt = mathRelOpComp2<FwNotEqual>(stack(), 2, s); break;
                        case Rpn::OpLess:
                            coreSt = mathRelOpComp2<FwLess>(stack(), 2, s); break;
                        case Rpn::OpLessEq:
                            coreSt = mathRelOpComp2<FwLessEq>(stack(), 2, s); break;
                        case Rpn::OpGreater:
                            coreSt = mathRelOpComp2<FwGreater>(stack(), 2, s); break;
                        case Rpn::OpGreaterEq:
                            coreSt = mathRelOpComp2<FwGreaterEq>(stack(), 2, s); break;
                        case Rpn::OpBitAnd:
                            coreSt = mathBitOpComp2<FwBitAnd>(stack(), 2, s); break;
                        case Rpn::OpBitOr:
                            coreSt = mathBitOpComp2<FwBitOr>(stack(), 2, s); break;
                        case Rpn::OpBitExor:
                            coreSt = mathBitOpComp2<FwBitExor>(stack(), 2, s); break;
                        case Rpn::OpBitShiftR:
                            coreSt = mathBitOpComp2<FwBitShiftRight>(stack(), 2, s); break;
                        case Rpn::OpBitShiftL:
                            coreSt = mathBitOpComp2<FwBitShiftLeft>(stack(), 2, s); break;
                        case Rpn::OpBitNot:
                            coreSt = mathBitOpComp1<FwBitNot>(stack(), 1, s); break;
                        case Rpn::OpNot:
                            coreSt = mathLogicOp1<FwNot>(stack(), 1, s); break;
                        case Rpn::OpAnd:
                        case Rpn::OpOr:
                            // Short circuit operations are never executed, 
                            // they are there only for expession formatting.
                            coreSt = true;
                            break;
                        case Rpn::OpSelect:
                            coreSt = mathFuncSelector2(stack(), 2, s); break;
                        default:
                            s.set( Status::Unsupported, "Unsupported operator."); 
                    }
                    if (!coreSt) {
                        appendLocation(s, loc);
                        return false;
                    }
                } catch (std::exception& exc) {
                    s.set(Status::Exception, std::string(exc.what())+" detected.");
                    appendLocation(s, loc);
                    return false;        
                }
                break;
            case Rpn::TFunctionCall: {
                // Find builtin
                auto& f = e->get<Rpn::FunctionCall>();
                auto* builtinPtr = contextStack_.getBuiltin(f.name);
                if (!builtinPtr) {
                    s.set(
                        Status::NotFound, 
                        std::string("Function '")+std::string(f.name)+"' not found."
                    );
                    appendLocation(s, loc);
                }
                // Check arity
                if (f.arity>builtinPtr->maxArity) {
                    s.set(
                        Status::BadArguments, 
                        std::string("Function '")+std::string(f.name)+"' accepts at most "+std::to_string(builtinPtr->maxArity)+" argument(s))."
                    );
                    appendLocation(s, loc);
                }
                if (f.arity<builtinPtr->minArity) {
                    s.set(
                        Status::BadArguments, 
                        std::string("Function '")+std::string(f.name)+"' requires at lest "+std::to_string(builtinPtr->maxArity)+" argument(s))."
                    );
                    appendLocation(s, loc);
                }
                // Call
                if (!builtinPtr->func(stack_, f.arity, s)) {
                    appendLocation(s, loc);
                    return false;
                }
                break;
            }
            case Rpn::TPackVec: {
                auto& p = e->get<Rpn::PackVec>();
                // Call
                if (!vectorPack(stack_, p.arity, s)) {
                    appendLocation(s, loc);
                    return false;
                }
                break;
            }
            case Rpn::TPackList: {
                auto& p = e->get<Rpn::PackList>();
                // Call
                if (!listPack(stack_, p.arity, s)) {
                    appendLocation(s, loc);
                    return false;
                }
                break;
            }
            case Rpn::TMergeList: {
                auto& p = e->get<Rpn::MergeList>();
                // Call
                if (!listMerge(stack_, p.arity, s)) {
                    appendLocation(s, loc);
                    return false;
                }
                break;
            }
            case Rpn::TJump: {
                auto& p = e->get<Rpn::Jump>();
                // Get offset, check for infinite loop
                if (p.offset==0) {
                    s.set(Status::Unsupported, "Internal error: infinite loop detected.");
                    appendLocation(s, loc);
                    return false;
                } 
                // Do not check destination, jump
                e += p.offset;
                isJump = true;
                break;
            }
            case Rpn::TBranch: {
                auto& br = e->get<Rpn::Branch>();
                bool negate = br.flags & Rpn::BrFalse; 
                bool keepOnBranch = br.flags&Rpn::BrKeepOnBranch;
                bool keepOnNoBranch = br.flags&Rpn::BrKeepOnNoBranch;
                auto offs = br.offset;
                // Get offset, check for infinite loop
                if (offs==0) {
                    s.set(Status::Unsupported, "Internal error: infinite loop detected.");
                    appendLocation(s, loc);
                    return false;
                } 

                // Get condition
                auto [ok, cond] = checkCondition();
                if (!ok) {
                    s.set(Status::Unsupported, "Don't know how to interpret value as a boolean.");
                    appendLocation(s, loc);
                    return false;
                }
                
                // Do not check destination, jump
                auto branch = cond^negate;
                if (cond^negate) {
                    e += offs;
                    isJump = true;
                }

                // Pop condition
                if (!(branch && keepOnBranch || !branch && keepOnNoBranch)) {
                    stack_.pop();
                }
                break;
            }
            case Rpn::TMakeBoolean: {
                auto [ok, cond] = checkCondition();
                if (!ok) {
                    s.set(Status::Unsupported, "Don't know how to interpret value as a boolean.");
                    appendLocation(s, loc);
                    return false;
                }
                
                // Pop value
                stack_.pop();

                // Replace with 1/0
                stack_.push(Value(cond?1:0));
                break;
            }
            default:
                s.set(Status::Unsupported, "RPN entry not supported.");
                appendLocation(s, loc);
                return false;
        }
        
        if (!isJump) {
            ++e;
        }
    }
    // Check stack depth
    if (stack_.size()==1) {
        stack_.pop(result);
    } else {
        DBGCHECK(true, std::string("Internal error. Final RPN stack contains ")+std::to_string(stack_.size())+" values insted of 1.");
        return false;
    }
    return true;
}

}
