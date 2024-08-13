#ifndef __FLAGS_DEFINED
#define __FLAGS_DEFINED

#include <type_traits>
#include "common.h"

namespace NAMESPACE {

// Define non-template flag operators
#define DEFINE_FLAG_OPERATORS(EnumType) \
    constexpr EnumType operator~(const EnumType a) { \
        return static_cast<EnumType>( \
            ~static_cast<std::underlying_type_t<EnumType>>(a) \
        ); \
    } \
    \
    constexpr EnumType operator|(const EnumType a, const EnumType b) { \
        return static_cast<EnumType>( \
            static_cast<std::underlying_type_t<EnumType>>(a) | \
            static_cast<std::underlying_type_t<EnumType>>(b) \
        ); \
    } \
    \
    constexpr EnumType operator&(const EnumType a, const EnumType b) { \
        return static_cast<EnumType>( \
            static_cast<std::underlying_type_t<EnumType>>(a) & \
            static_cast<std::underlying_type_t<EnumType>>(b) \
        ); \
    } \
    constexpr EnumType operator<<(const EnumType a, const int b) { \
        return static_cast<EnumType>( \
            static_cast<std::underlying_type_t<EnumType>>(a) << b \
        ); \
    } \
    constexpr EnumType operator>>(const EnumType a, const int b) { \
        return static_cast<EnumType>( \
            static_cast<std::underlying_type_t<EnumType>>(a) >> b \
        ); \
    } \
    constexpr bool operator==(const EnumType a, const int b) { \
        return static_cast<std::underlying_type_t<EnumType>>(a) == b; \
    } \
    constexpr bool operator!=(const EnumType a, const int b) { \
        return static_cast<std::underlying_type_t<EnumType>>(a) != b; \
    } 

template<typename EnumType> class FlagBase {
public:
    using Flags = EnumType;

    const static EnumType NoFlags = EnumType(0);
    
    FlagBase() : flags_(EnumType(0)) {};
    FlagBase(EnumType f) : flags_(f) {};

    constexpr EnumType flags() { return flags_; };
    constexpr EnumType maskedFlags(Flags mask) { 
        return static_cast<EnumType>(
            static_cast<NumType>(flags_) & static_cast<NumType>(mask)); 
    };
    constexpr bool checkFlags(Flags f) { return maskedFlags(f) == f; }; 
    void setFlags(Flags f) { 
        flags_ = static_cast<EnumType>(static_cast<NumType>(flags_) | static_cast<NumType>(f)); 
    };
    void clearFlags(Flags f = static_cast<EnumType>(~NumType(0))) {
        flags_ = static_cast<EnumType>(static_cast<NumType>(flags_) & ~static_cast<NumType>(f)); 
    };
    
protected:
    using NumType = typename std::underlying_type<EnumType>::type;
    EnumType flags_;
};

}

#endif
