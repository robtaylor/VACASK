#ifndef __RPNEXPR_DEFINED
#define __RPNEXPR_DEFINED

#include <variant>
#include <limits>
#include "value.h"
#include "identifier.h"
#include "common.h"


namespace NAMESPACE {

class Rpn {
public:
    enum OpCode : char { 
        OpPlus, OpMinus, OpTimes, OpDivide, OpUMinus, OpPower, 
        OpEqual, OpNotEqual, OpLess, OpLessEq, OpGreater, OpGreaterEq, 
        OpBitAnd, OpBitOr, OpBitExor, OpBitNot, OpBitShiftR, OpBitShiftL, 
        OpAnd, OpOr, OpNot, OpQuestion, 
        OpSelect
    };
    enum BrFlags : char {
        BrKeepOnBranch=1, BrKeepOnNoBranch=2, BrFalse=4, BrHidden=8
    };

    static const LocationIndex badLocationIndex = std::numeric_limits<LocationIndex>::max();
    static const LocationIndex maxLocationIndex = std::numeric_limits<LocationIndex>::max() - 1;
    typedef RpnArity Arity;
    typedef RpnJumpOffset JumpOffset;
    // 1+8 = 9 bytes
    typedef struct Op {
        OpCode code;
        LocationIndex loc;
        Op(const OpCode& c) : code(c), loc(badLocationIndex) {};
        bool operator==(OpCode op) { return code==op; };
    } Op;
    // 4+8 = 12 bytes
    typedef struct Identifier {
        Id name;
        LocationIndex loc;
        Identifier(const std::string&& s) : name(std::move(s)), loc(badLocationIndex) {};
    } Identifier;
    // 4+4+8 = 16 bytes
    typedef struct FunctionCall {
        Id name;
        Arity arity;
        LocationIndex loc;
        FunctionCall(const std::string&& s, Arity a) : name(std::move(s)), arity(a), loc(badLocationIndex) {};
    } FunctionCall;
    // 4+4 = 8 bytes
    typedef struct PackVec {
        Arity arity;
        LocationIndex loc;
        PackVec(Arity a) : arity(a), loc(badLocationIndex) {};
    } PackVec;
    // 4+4 = 8 bytes
    typedef struct PackList {
        Arity arity;
        LocationIndex loc;
        PackList(Arity a) : arity(a), loc(badLocationIndex) {};
    } PackList;
    // 4+4 = 8 bytes
    typedef struct MergeList {
        Arity arity;
        LocationIndex loc;
        MergeList(Arity a) : arity(a), loc(badLocationIndex) {};
    } MergeList;
    // 4+8 = 12 bytes
    typedef struct Jump {
        JumpOffset offset;
        BrFlags flags;
        LocationIndex loc;
        Jump(JumpOffset o, BrFlags f) : offset(o), flags(f), loc(badLocationIndex) {};
    } Jump;
    // 4+8 = 12 bytes
    typedef struct Branch {
        JumpOffset offset;
        BrFlags flags;
        LocationIndex loc;
        Branch(JumpOffset o, BrFlags f) : offset(o), flags(f), loc(badLocationIndex) {};
    } Branch;
    // 8 = 8 bytes
    typedef struct MakeBoolean {
        LocationIndex loc;
        MakeBoolean() : loc(badLocationIndex) {};
    } MakeBoolean;
    
    enum Type : char { 
        TValue=0, TOp=1, TIdentifier=2, TFunctionCall=3, 
        TPackVec=4, TPackList=5, TMergeList=6, 
        TJump=7, TBranch=8, TMakeBoolean=9
    };

    class Entry {
    public:
        Entry           (const Entry&)  = delete;
        Entry           (      Entry&&) = default;
        Entry& operator=(const Entry&)  = delete;
        Entry& operator=(      Entry&&) = default;

        template<typename T> Entry(T&& other) : data(std::move(other)) {};
        
        Type type() const { return Type(data.index()); };
        void setLocation(LocationIndex li) {
            switch (type()) {
                case TOp: std::get<Op>(data).loc = li; break;
                case TIdentifier: std::get<Identifier>(data).loc = li; break;
                case TFunctionCall: std::get<FunctionCall>(data).loc = li; break;
                case TPackVec: std::get<PackVec>(data).loc = li; break;
                case TPackList: std::get<PackList>(data).loc = li; break;
                case TMergeList: std::get<MergeList>(data).loc = li; break;
                case TJump: std::get<Jump>(data).loc = li; break;
                case TBranch: std::get<Branch>(data).loc = li; break;
                case TMakeBoolean: std::get<MakeBoolean>(data).loc = li; break;
            }
        };
        LocationIndex location() const {
            switch (type()) {
                case TOp: return std::get<Op>(data).loc;
                case TIdentifier: return std::get<Identifier>(data).loc;
                case TFunctionCall: return std::get<FunctionCall>(data).loc;
                case TPackVec: return std::get<PackVec>(data).loc;
                case TPackList: return std::get<PackList>(data).loc;
                case TMergeList: return std::get<MergeList>(data).loc;
                case TJump: return std::get<Jump>(data).loc;
                case TBranch: return std::get<Branch>(data).loc;
                case TMakeBoolean: return std::get<MakeBoolean>(data).loc;
                default: return badLocationIndex;
            }
        };
        template<typename T> T& get() { return std::get<T>(data); };
        template<typename T> const T& get() const { return std::get<T>(data); };
        
        std::variant<Value, Op, Identifier, FunctionCall, PackVec, PackList, MergeList, Jump, Branch, MakeBoolean> data;
    };
    
    typedef std::vector<Entry> Expression;
    typedef std::vector<Loc> Locations;
    
    Rpn();

    Rpn           (const Rpn&)  = delete;
    Rpn           (      Rpn&&) = default;
    Rpn& operator=(const Rpn&)  = delete;
    Rpn& operator=(      Rpn&&) = default;
    
    inline auto begin() { return expr.begin(); };
    inline auto end() { return expr.end(); };
    inline auto cbegin() const { return expr.cbegin(); };
    inline auto cend() const { return expr.cend(); };

    inline bool endsWithMakeBoolean() {
        return expr.size()!=0 && expr.back().type()==Rpn::TMakeBoolean;
    };
    inline void extend(Rpn &&other) {
        auto locBase = locations.size();
        for(auto it=other.begin(); it!=other.end(); ++it) {
            it->setLocation(it->location()+locBase);
            expr.push_back(std::move(*it));
        }
        if (locBase+other.locations.size()>maxLocationIndex) {
            throw std::length_error("Too many location descriptors. Expression is too long."); 
        }
        for(auto it=other.locations.begin(); it!=other.locations.end(); ++it) {
            locations.push_back(std::move(*it));
        }
    };
    inline void extend(Entry&& other, Loc l) { 
        auto locIndex = locations.size();
        if (locIndex>maxLocationIndex) {
            throw std::length_error("Too many location descriptors. Expression is too long."); 
        }
        locations.push_back(l);
        other.setLocation(locIndex);
        expr.push_back(std::move(other)); 
    };
    inline const Loc& location(const Entry& e) const { 
        auto li = e.location();
        if (li!=badLocationIndex) {
            return locations[li];
        } else {
            return Loc::bad;
        }
    }
    
    inline size_t size() const noexcept { return expr.size(); }

    inline const Rpn::Entry& operator[](size_t i) const { return expr[i]; }

    std::string str() const;
    
    friend std::ostream& operator<<(std::ostream& os, const Rpn& expr);

private:
    Expression expr;
    Locations locations;
    static std::unordered_map<OpCode, std::tuple<const char*, int>> opMap;

    std::string parenthesize(const std::string& s) const {
        return std::string("(")+s+")";
    };
};

DEFINE_FLAG_OPERATORS(Rpn::BrFlags);

}

#endif
