#include "rpnexpr.h"
#include "common.h"


namespace NAMESPACE {

Rpn::Rpn() {
}

// This replicates the operator priority from the parser
// Higher number = higher priority
std::unordered_map<Rpn::OpCode, std::tuple<const char*, int>> Rpn::opMap = {
    { OpOr, {"||", -7} },  
    { OpAnd, {"&&", -6} },  
    { OpBitOr, {"|", -5} },  
    { OpBitExor, {"^", -4} },  
    { OpBitAnd, {"&", -3} },  
    { OpEqual, {"==", -2} }, { OpNotEqual, {"!=", -2} }, 
    { OpLess, {"<", -1} }, { OpLessEq, {"<=", -1} }, { OpGreater, {">", -1} },  { OpGreaterEq, {">=", -1} }, 
    { OpBitShiftL, {"<<", 0 } },  { OpBitShiftR, {">>", 0} }, 
    { OpPlus, {"+", 1} }, { OpMinus, {"-", 1} }, 
    { OpTimes, {"*", 2} }, { OpDivide, {"/", 2} }, 
    { OpPower, {"**", 3} }, 
    { OpUMinus, {"-", 4} }, { OpNot, {"!", 4} }, { OpBitNot, {"~", 4} }, 
    { OpSelect, {"[]", 5} }, 
};

std::string Rpn::str() const {
    std::vector<std::tuple<std::string, int>> sstack;
    OpCode code;
    for(auto e=expr.cbegin(); e!=expr.cend(); ++e) {
        switch (e->type()) {
            case Rpn::TValue: {
                const Value& v = e->get<Value>();
                sstack.push_back({v.str(), 1000});
                continue;
            }
            case Rpn::TIdentifier: {
                sstack.push_back({std::string(e->get<Identifier>().name), 1000});
                continue;
            }
            case Rpn::TOp: {
                code = e->get<Op>().code;
                auto it = opMap.find(code);
                int priority = -1000;
                const char* opStr = "";
                if (it!=opMap.end()) {
                    std::tie(opStr, priority) = it->second;
                }
                switch (code) {
                    case OpUMinus:
                    case OpNot:
                    case OpBitNot: {
                        // Unary prefix operators
                        auto [ ex, exprPr ] = std::move(sstack.back());
                        sstack.pop_back();
                        if (exprPr<priority) {
                            // expression has lower prority, use -(ex)
                            sstack.push_back({std::move(std::string(opStr)+"("+ex+")"), exprPr});
                        } else {
                            // expression has higher prority, use -ex
                            sstack.push_back({std::move(std::string(opStr)+ex), exprPr});
                        }
                        break;
                    }
                    case OpPlus: 
                    case OpMinus: 
                    case OpTimes: 
                    case OpDivide: 
                    case OpPower: 
                    case OpBitShiftL: 
                    case OpBitShiftR: 
                    case OpLess: 
                    case OpLessEq: 
                    case OpGreater: 
                    case OpGreaterEq: 
                    case OpEqual: 
                    case OpNotEqual: 
                    case OpBitAnd: 
                    case OpBitOr: 
                    case OpBitExor: 
                    case OpAnd: 
                    case OpOr: {
                        // Binary infix operators 
                        std::string ex1, ex2;
                        int exprPr1, exprPr2;
                        std::tie(ex2, exprPr2) = std::move(sstack.back());
                        sstack.pop_back();
                        std::tie(ex1, exprPr1) = std::move(sstack.back());
                        sstack.pop_back();
                        sstack.push_back({
                            (exprPr1<priority ? parenthesize(ex1) : ex1) +
                            opStr + 
                            (exprPr2<priority ? parenthesize(ex2) : ex2), 
                            priority
                        });
                        break;
                    }
                    case OpSelect: {
                        // Selector
                        std::string ex1, ex2;
                        int exprPr1, exprPr2;
                        std::tie(ex1, exprPr1) = std::move(sstack.back());
                        sstack.pop_back();
                        std::tie(ex2, exprPr2) = std::move(sstack.back());
                        sstack.pop_back();
                        sstack.push_back({
                            (exprPr1<priority ? parenthesize(ex1) : ex1) +
                            "["+ex2+"]", 
                            priority
                        });
                        break;
                    }
                }
                break;
            }
            case TFunctionCall: {
                // f()
                // f(x1)
                // f(x1, x2, ..., xn)
                auto name = e->get<FunctionCall>().name;
                auto n = e->get<FunctionCall>().arity;
                std::string txt = std::string(name)+"(";
                auto j = sstack.size()-n;
                for(decltype(n) i=0; i<n; i++, j++) {
                    if (i>0 && n>1)  {
                        txt += ", ";
                    }
                    txt += std::get<0>(sstack.at(j));
                }
                txt+=")"; 
                sstack.resize(sstack.size()-n);
                sstack.push_back({std::move(txt), 1000});
                break;
            }
            case TPackVec: {
                // [], [,]
                // [x1]
                // [x1, x2, ..., xn]
                auto n = e->get<PackVec>().arity;
                std::string txt = "[";
                auto j = sstack.size()-n;
                for(decltype(n) i=0; i<n; i++, j++) {
                    if (i>0 && n>1)  {
                        txt += ", ";
                    }
                    txt += std::get<0>(sstack.at(j));
                }
                txt+="]"; 
                sstack.resize(sstack.size()-n);
                sstack.push_back({std::move(txt), 1000});
                break;
            }
            case TPackList: {
                // [;]
                // [x1;]
                // [x1; x2; ...; xn]
                auto n = e->get<PackVec>().arity;
                std::string txt = "[";
                auto j = sstack.size()-n;
                for(decltype(n) i=0; i<n; i++, j++) {
                    if (i>0 && n>1)  {
                        txt += "; ";
                    }
                    txt += std::get<0>(sstack.at(j));
                }
                if (n<2) {
                    txt +=";";
                }
                txt+="]"; 
                sstack.resize(sstack.size()-n);
                sstack.push_back({std::move(txt), 1000});
                break;
            }
            case TMergeList: {
                // [:]
                // [x1:]
                // [x1: x2: ...: xn]
                auto n = e->get<MergeList>().arity;
                std::string txt = "[";
                auto j = sstack.size()-n;
                for(decltype(n) i=0; i<n; i++, j++) {
                    if (i>0 && n>1)  {
                        txt += ": ";
                    }
                    txt += std::get<0>(sstack.at(j));
                }
                if (n<2) {
                    txt +=":";
                }
                txt+="]"; 
                sstack.resize(sstack.size()-n);
                sstack.push_back({std::move(txt), 1000});
                break;
            }
            case TMakeBoolean:
            case TJump: 
            case TBranch: 
                // Do not format if hidden
                break;
            default:
                return "<expr UNSUPPORTED>";
        }
    }
    return std::get<0>(sstack.back()); 
}

std::ostream& operator<<(std::ostream& os, const Rpn& expr) {
    os << expr.str();
    return os;
}

}
