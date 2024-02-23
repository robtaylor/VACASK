#ifndef __VALUE_DEFINED
#define __VALUE_DEFINED

#include <iostream>
#include <string>
#include <vector>
#include <variant>
#include <type_traits>
#include "filestack.h"
#include "sourceloc.h"
#include "status.h"
#include "common.h"
#include "flags.h"
#include <cinttypes>


namespace NAMESPACE {

class Value;

// Can set Int to int, int32_t, or int64_t
typedef IntegerValue Int;
typedef double Real;
typedef std::string String;
typedef std::vector<Int> IntVector;
typedef std::vector<Real> RealVector;
typedef std::vector<String> StringVector;
typedef std::vector<Value> ValueVector; // List of values

enum class ValueType : char { 
    Zero=0, VectorBit=8, 
    Int=0, Real=1, String=2, 
    Value=4,  // Never used in Value class for storage
              // Used only in introspection mechanism
    IntVec=8, RealVec=9, StringVec=10, ValueVec=12 
};
DEFINE_FLAG_OPERATORS(ValueType);

class Value {
public:
    // Value and Value vector are used unly for access to introspected structures. 
    // At this point variant and variant vector are not an option for Value type. 
    // In the furture variant vector might be implemented. 
    using Type = ValueType;
    
    // Value           (const Value&)  = delete;
    // Value           (      Value&&) = default;
    // Value& operator=(const Value&)  = delete;
    // Value& operator=(      Value&&) = default;
    
    // Copy constructor
    inline Value(const Value& other) {
        switch (other.type_) {
            case Type::Int:        iVal = other.iVal; break;
            case Type::Real:       rVal = other.rVal; break;
            case Type::String:     sVal = new String(*(other.sVal)); break;
            case Type::IntVec:     iVec = new IntVector(*(other.iVec)); break;
            case Type::RealVec:    rVec = new RealVector(*(other.rVec)); break;
            case Type::StringVec:  sVec = new StringVector(*(other.sVec)); break;
            case Type::ValueVec: vVec = new ValueVector(*(other.vVec)); break;
        }
        type_ = other.type_;
    };

    // Move constructor
    inline Value(Value &&other) noexcept {
        switch (other.type_) {
            case Type::Int:        iVal = other.iVal; break;
            case Type::Real:       rVal = other.rVal; break;
            case Type::String:     sVal = other.sVal; other.sVal = nullptr; break;
            case Type::IntVec:     iVec = other.iVec; other.iVec = nullptr; break;
            case Type::RealVec:    rVec = other.rVec; other.rVec = nullptr; break;
            case Type::StringVec:  sVec = other.sVec; other.sVec = nullptr; break;
            case Type::ValueVec: vVec = other.vVec; other.vVec = nullptr; break;
        }  
        type_ = other.type_;
    }

    // Constructors from values passed by value (copy)
    inline Value(const Int &other) : iVal(other), type_(Type::Int) {};
    inline Value(const Real &other) : rVal(other), type_(Type::Real) {};
    inline Value(const String &other) : sVal(new String(other)), type_(Type::String) {};
    inline Value(const IntVector &other) : iVec(new IntVector(other)), type_(Type::IntVec) {};
    inline Value(const RealVector &other) : rVec(new RealVector(other)), type_(Type::RealVec) {};
    inline Value(const StringVector &other) : sVec(new StringVector(other)), type_(Type::StringVec) {};
    inline Value(const ValueVector &other) : vVec(new ValueVector(other)), type_(Type::ValueVec) {};

    // Constructors from values passed by rvalue reference (move)
    inline Value(String &&other) noexcept : type_(Type::String) { sVal = new String(std::move(other)); };
    inline Value(IntVector &&other) noexcept : type_(Type::IntVec) { iVec = new IntVector(std::move(other)); };
    inline Value(RealVector &&other) noexcept : type_(Type::RealVec) { rVec = new RealVector(std::move(other)); };
    inline Value(StringVector &&other) noexcept : type_(Type::StringVec) { sVec = new StringVector(std::move(other)); };
    inline Value(ValueVector &&other) noexcept : type_(Type::ValueVec) { vVec = new ValueVector(std::move(other)); };

    // Default constructor, creates zero Int
    inline Value() noexcept : Value(Int(0)) {};

    // Destroy value
    inline ~Value() { 
        switch (type_) {
            case Type::String: delete sVal; break;
            case Type::IntVec: delete iVec; break;
            case Type::RealVec: delete rVec; break;
            case Type::StringVec: delete sVec; break;
            case Type::ValueVec: delete vVec; break;
        }
    }

    // Extract value, returns reference to value
    template<typename T> T& val();
    template<typename T> T& val() const;
    
    // Convert
    bool convert(Type to, Value& dest, Status& s=Status::ignore);
    bool convertInPlace(Type to, Status& s=Status::ignore);

    // Size of vector, return -1 on non-vectors
    // Nonnegative for indexable values
    inline size_t size() const noexcept {
        switch (type_) {
            case Type::IntVec: return iVec->size();
            case Type::RealVec: return rVec->size();
            case Type::StringVec: return sVec->size();
            case Type::ValueVec: return vVec->size();
            default: return -1;
        }
        return -1; // Reaching this is an error
    };

    // Check type
    inline Type type() const noexcept { return type_; };
    inline Type scalarType() const noexcept { return Type(type_ & ~Type::VectorBit); }; 
    inline bool isVector() const noexcept { return (type_ & Type::VectorBit) != 0; };

    // Swap value with other
    friend inline void swap(Value &a, Value &b) {
        if (&a == &b) {
            return;
        }
        Value temp = std::move(a);
        a.~Value(); new(&a) Value(std::move(b));
        b.~Value(); new(&b) Value(std::move(temp));
    };

    // Copy assignment
    inline Value &operator=(const Value& other) {
        Value temp = other; // Make a copy
        swap(*this, temp);
        return *this;
    };
    inline Value &operator=(Int other) { this->~Value(); iVal = other; type_ = Type::Int; return *this; };
    inline Value &operator=(Real other) { this->~Value(); rVal = other; type_ = Type::Real; return *this; };
    inline Value &operator=(const String& other) { this->~Value(); sVal = new String(other); type_ = Type::String; return *this; };
    inline Value &operator=(const IntVector& other) { this->~Value(); iVec = new IntVector(other); type_ = Type::IntVec; return *this; };
    inline Value &operator=(const RealVector& other) { this->~Value(); rVec = new RealVector(other); type_ = Type::RealVec; return *this; };
    inline Value &operator=(const StringVector& other) { this->~Value(); sVec = new StringVector(other); type_ = Type::StringVec; return *this; };
    inline Value &operator=(const ValueVector& other) { this->~Value(); vVec = new ValueVector(other); type_ = Type::ValueVec; return *this; };

    // Move assignment
    inline Value &operator=(Value &&other) noexcept {
        this->~Value(); new(this) Value(std::move(other));
        return *this;
    };
    inline Value &operator=(String&& other) { 
        if (type_!=Type::String) { this->~Value(); sVal = new String(std::move(other)); type_ = Type::String; } 
        else { *sVal = std::move(other); }
        return *this;
    };
    inline Value &operator=(IntVector&& other) { 
        if (type_!=Type::IntVec) { this->~Value(); iVec = new IntVector(std::move(other)); type_ = Type::IntVec; } 
        else { *iVec = std::move(other); }
        return *this;
    };
    inline Value &operator=(RealVector&& other) { 
        if (type_!=Type::RealVec) { this->~Value(); rVec = new RealVector(std::move(other)); type_ = Type::RealVec; } 
        else { *rVec = std::move(other); }
        return *this;
    };
    inline Value &operator=(StringVector&& other) { 
        if (type_!=Type::StringVec) { this->~Value(); sVec = new StringVector(std::move(other)); type_ = Type::StringVec; } 
        else { *sVec = std::move(other); }
        return *this;
    };
    inline Value &operator=(ValueVector&& other) { 
        if (type_!=Type::ValueVec) { this->~Value(); vVec = new ValueVector(std::move(other)); type_ = Type::ValueVec; } 
        else { *vVec = std::move(other); }
        return *this;
    };

    // Extract scalar from vector
    bool getScalar(Value& v, size_t ndx, Status& s=Status::ignore) const;

    // Compare
    friend bool operator==(const Value &, const Value &);
    friend bool operator!=(const Value &, const Value &);
    
    // Output 
    friend std::ostream& operator<<(std::ostream& os, const Value& obj);
    std::string str() const;

    // Instance type name
    std::string typeName() { return typeCodeToName(type_); };

    // Type code from type
    template<typename T> static Type typeCode();

    // Type code to name
    static std::string typeCodeToName(Type t);

    // typeCode to type
    // Define a templated struct and specialize it for all available types
    template<Type I> struct TypeTag_ {};
    // Now, define a new type
    template<Type I> using CodeToType = typename TypeTag_<I>::type;

    // Type to scalar type
    // Scalar type of Value can be anything, should cause compile error
    template<typename T> struct ScalarType_ { 
        static_assert(std::is_same<T, Value>::value, "Cannot determine scalar type for Value."); 
    };
    template<typename T> using ScalarType = typename ScalarType_<T>::type;

    // Is type a vector type
    // Value can be either a scalar or a vector, should cause compile error
    template<typename T> struct IsVectorType {
        static_assert(std::is_same<T, Value>::value, "Cannot determine if Value type is a vector."); 
    };
    
    // Type to vector type, default is for vectorize = false
    // Specializations given for vectorize = true
    template<typename T, bool vectorize> struct VectorType_ { using type = T; };
    template<typename T, bool vectorize> using VectorType = typename VectorType_<T,vectorize>::type;

private:
    union {
        Int iVal;
        Real rVal;
        String* sVal;
        IntVector* iVec;
        RealVector* rVec;
        StringVector* sVec;
        ValueVector* vVec;
    };
    Type type_;
};

// Helper function that is computed at compile time
constexpr uint32_t type_pair(Value::Type t1, Value::Type t2) {
    return (static_cast<uint32_t>(t1)<<16)+static_cast<uint32_t>(t2);
}

// Value reference retrievers
template<> inline Int& Value::val<Int>() { return iVal; }
template<> inline Real& Value::val<Real>() { return rVal; }
template<> inline String& Value::val<String>() { return *sVal; }
template<> inline IntVector& Value::val<IntVector>() { return *iVec; }
template<> inline RealVector& Value::val<RealVector>() { return *rVec; }
template<> inline StringVector& Value::val<StringVector>() { return *sVec; }
template<> inline ValueVector& Value::val<ValueVector>() { return *vVec; }

template<> inline const Int& Value::val<const Int>() const { return iVal; }
template<> inline const Real& Value::val<const Real>() const { return rVal; }
template<> inline const String& Value::val<const String>() const { return *sVal; }
template<> inline const IntVector& Value::val<const IntVector>() const { return *iVec; }
template<> inline const RealVector& Value::val<const RealVector>() const { return *rVec; }
template<> inline const StringVector& Value::val<const StringVector>() const { return *sVec; }
template<> inline const ValueVector& Value::val<const ValueVector>() const { return *vVec; }

// Conversion of type to type code
template<> inline Value::Type Value::typeCode<Int>() { return Value::Type::Int; }
template<> inline Value::Type Value::typeCode<Real>() { return Value::Type::Real; }
template<> inline Value::Type Value::typeCode<String>() { return Value::Type::String; }
template<> inline Value::Type Value::typeCode<IntVector>() { return Value::Type::IntVec; }
template<> inline Value::Type Value::typeCode<RealVector>() { return Value::Type::RealVec; }
template<> inline Value::Type Value::typeCode<StringVector>() { return Value::Type::StringVec; }

template<> inline Value::Type Value::typeCode<Value>() { return Value::Type::Value; }
template<> inline Value::Type Value::typeCode<ValueVector>() { return Value::Type::ValueVec; }

// Conversion of type code to type
// Specialize Value::TypeTag_ for all available types
template<> struct Value::TypeTag_<Value::Type::Int> { using type=Int; };
template<> struct Value::TypeTag_<Value::Type::Real> { using type=Real; };
template<> struct Value::TypeTag_<Value::Type::String> { using type=String; };
template<> struct Value::TypeTag_<Value::Type::IntVec> { using type=IntVector; };
template<> struct Value::TypeTag_<Value::Type::RealVec> { using type=RealVector; };
template<> struct Value::TypeTag_<Value::Type::StringVec> { using type=StringVector; };
template<> struct Value::TypeTag_<Value::Type::ValueVec> { using type=ValueVector; }; 

// Conversion of type to scalar type
// Specialize Value::ScalarType for all available types
template<> struct Value::ScalarType_<Int> { using type=Int; };
template<> struct Value::ScalarType_<Real> { using type=Real; };
template<> struct Value::ScalarType_<String> { using type=String; };
template<> struct Value::ScalarType_<IntVector> { using type=Int; };
template<> struct Value::ScalarType_<RealVector> { using type=Real; };
template<> struct Value::ScalarType_<StringVector> { using type=String; };
template<> struct Value::ScalarType_<ValueVector> { using type=Value; };

// Convert type to a vector type
template<> struct Value::VectorType_<Int,true> { using type=IntVector; };
template<> struct Value::VectorType_<Real,true> { using type=RealVector; };
template<> struct Value::VectorType_<String,true> { using type=StringVector; };
template<> struct Value::VectorType_<Value,true> { using type=ValueVector; };

// Check for vector type
template<> struct Value::IsVectorType<Int> { static const bool value=false; };
template<> struct Value::IsVectorType<Real> { static const bool value=false; };
template<> struct Value::IsVectorType<String> { static const bool value=false; };
template<> struct Value::IsVectorType<IntVector> { static const bool value=true; };
template<> struct Value::IsVectorType<RealVector> { static const bool value=true; };
template<> struct Value::IsVectorType<StringVector> { static const bool value=true; };
template<> struct Value::IsVectorType<ValueVector> { static const bool value=true; };

}

#endif
