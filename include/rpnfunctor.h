#ifndef __RPNFUNCTORS_DEFINED
#define __RPNFUNCTORS_DEFINED

#include "value.h"
#include "status.h"
#include <cmath>
#include <type_traits>
#include "libplatform.h"
#include "common.h"


namespace NAMESPACE {

// Functors for wrapping mathematical functions so that we can pass them to templates
// and the compiler can then inline them. 
// If we pass a simple function pointer the compiler does not know how to inline it. 

// All functors have an ok() member that tests the validity of arguments

struct FwUminus { 
    template<typename T> T operator()(T x) { return -x; }; 
    template<typename T> bool ok(T x, Status& s) { return true; }; 
};

struct FwPlus { 
    template<typename T1, typename T2> auto operator()(T1 x1, T2 x2) -> decltype(x1+x2) { return x1+x2; }; 
    template<typename T1, typename T2> bool ok(T1 x1, T2 x2, Status& s) { return true; }; 
};

struct FwMinus { 
    template<typename T1, typename T2> auto operator()(T1 x1, T2 x2) -> decltype(x1-x2) { return x1-x2; }; 
    template<typename T1, typename T2> bool ok(T1 x1, T2 x2, Status& s) { return true; }; 
};

struct FwTimes { 
    template<typename T1, typename T2> inline auto operator()(T1 x1, T2 x2) -> decltype(x1*x2) { return x1*x2; }; 
    template<typename T1, typename T2> inline bool ok(T1 x1, T2 x2, Status& s) { return true; }; 
};

struct FwDivide { 
    template<typename T1, typename T2> auto operator()(T1 x1, T2 x2) -> decltype(x1/x2) { return x1/x2; }; 
    template<typename T1, typename T2> bool ok(T1 x1, T2 x2, Status& s) { 
        if (x2!=0) {
            return true;
        } else {
            s.set(Status::DivZero, "Division by zero.");
            return false;
        }
    }; 
};

struct FwPower { 
    template<typename T1, typename T2> auto operator()(T1 x1, T2 x2) -> decltype(std::pow(x1, x2)) { return std::pow(x1, x2); }; 
    template<typename T1, typename T2> bool ok(T1 x1, T2 x2, Status& s) { 
        // TODO: check all bad cases of std::pow()
        if (x1==0 && x2<=0) {
            s.set(Status::DivZero, "Zero raised to negative or zero power.");
            return false;
        } else {
            return true;
        }
    }; 
};

struct FwEqual { 
    template<typename T1, typename T2> Int operator()(T1 x1, T2 x2) { return x1==x2; }; 
    template<typename T1, typename T2> bool ok(T1 x1, T2 x2, Status& s) { return true; }; 
};

struct FwNotEqual { 
    template<typename T1, typename T2> Int operator()(T1 x1, T2 x2) { return x1!=x2; }; 
    template<typename T1, typename T2> bool ok(T1 x1, T2 x2, Status& s) { return true; }; 
};

struct FwLess { 
    template<typename T1, typename T2> Int operator()(T1 x1, T2 x2) { return x1<x2; }; 
    template<typename T1, typename T2> bool ok(T1 x1, T2 x2, Status& s) { return true; }; 
};

struct FwLessEq { 
    template<typename T1, typename T2> Int operator()(T1 x1, T2 x2) { return x1<=x2; }; 
    template<typename T1, typename T2> bool ok(T1 x1, T2 x2, Status& s) { return true; }; 
};

struct FwGreater { 
    template<typename T1, typename T2> Int operator()(T1 x1, T2 x2) { return x1>x2; }; 
    template<typename T1, typename T2> bool ok(T1 x1, T2 x2, Status& s) { return true; }; 
};

struct FwGreaterEq { 
    template<typename T1, typename T2> Int operator()(T1 x1, T2 x2) { return x1>=x2; }; 
    template<typename T1, typename T2> bool ok(T1 x1, T2 x2, Status& s) { return true; }; 
};

struct FwBitAnd { 
    Int operator()(Int x1, Int x2) { return x1&x2; }; 
    bool ok(Int x1, Int x2, Status& s) { return true; }; 
};

struct FwBitOr { 
    Int operator()(Int x1, Int x2) { return x1|x2; }; 
    bool ok(Int x1, Int x2, Status& s) { return true; }; 
};

struct FwBitExor { 
    Int operator()(Int x1, Int x2) { return x1^x2; }; 
    bool ok(Int x1, Int x2, Status& s) { return true; }; 
};

struct FwBitNot{ 
    Int operator()(Int x) { return ~x; }; 
    bool ok(Int x, Status& s) { return true; }; 
};

struct FwBitShiftRight{ 
    Int operator()(Int x, Int n) { return x>>n; }; 
    bool ok(Int x, Int n, Status& s) { return true; }; 
};

struct FwBitShiftLeft{ 
    Int operator()(Int x, Int n) { return x<<n; }; 
    bool ok(Int x, Int n, Status& s) { return true; }; 
};


template<typename T> inline Int makeBool(T x) { return x ? 1 : 0; }; 
inline Int makeBool(const String& x) { return x.size() ? 1 : 0; }; 

struct FwAnd { 
    template<typename T1, typename T2> Int operator()(T1 x1, T2 x2) { return makeBool(x1)&&makeBool(x2); }; 
    template<typename T1, typename T2> bool ok(T1 x1, T2 x2, Status& s) { return true; }; 
};

struct FwOr { 
    template<typename T1, typename T2> Int operator()(T1 x1, T2 x2) { return makeBool(x1)||makeBool(x2); }; 
    template<typename T1, typename T2> bool ok(T1 x1, T2 x2, Status& s) { return true; }; 
};

struct FwNot { 
    template<typename T> Int operator()(T x) { return !makeBool(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true; }; 
};

struct FwSin { 
    template<typename T> auto operator()(T x) -> decltype(std::sin(x)) { return std::sin(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true;}; 
};

struct FwCos { 
    template<typename T> auto operator()(T x) -> decltype(std::cos(x)) { return std::cos(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true;}; 
};

struct FwTan { 
    template<typename T> auto operator()(T x) -> decltype(std::tan(x)) { return std::tan(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true;}; 
};
struct FwAsin { 
    template<typename T> auto operator()(T x) -> decltype(std::asin(x)) { return std::asin(x); }; 
    template<typename T> bool ok(T x, Status& s) { 
        if ((x>=-1) && (x<=1)) {
            return true;
        } else {
            s.set(Status::Domain, "Domain error.");
            return false;
        } 
    }; 
};

struct FwAcos { 
    template<typename T> auto operator()(T x) -> decltype(std::acos(x)) { return std::acos(x); }; 
    template<typename T> bool ok(T x, Status& s) { 
        if ((x>=-1) && (x<=1)) {
            return true;
        } else {
            s.set(Status::Domain, "Domain error.");
            return false;
        }
    }; 
};

struct FwAtan { 
    template<typename T> auto operator()(T x) -> decltype(std::atan(x)) { return std::atan(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true;}; 
};


struct FwSinh { 
    template<typename T> auto operator()(T x) -> decltype(std::sinh(x)) { return std::sinh(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true;}; 
};

struct FwCosh { 
    template<typename T> auto operator()(T x) -> decltype(std::cosh(x)) { return std::cosh(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true;}; 
};

struct FwTanh { 
    template<typename T> auto operator()(T x) -> decltype(std::tanh(x)) { return std::tanh(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true;}; 
};

struct FwAsinh { 
    template<typename T> auto operator()(T x) -> decltype(std::asinh(x)) { return std::asinh(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true;}; 
};

struct FwAcosh { 
    template<typename T> auto operator()(T x) -> decltype(std::acosh(x)) { return std::acosh(x); }; 
    template<typename T> bool ok(T x, Status& s) { 
        if (x>=1) {
            return true;
        } else {
            s.set(Status::Domain, "Domain error.");
            return false;
        }
    }; 
};

struct FwAtanh { 
    template<typename T> auto operator()(T x) -> decltype(std::atanh(x)) { return std::atanh(x); }; 
    template<typename T> bool ok(T x, Status& s) { 
        if ((x>=-1) && (x<=1)) {
            return true;
        } else {
            s.set(Status::Domain, "Domain error.");
            return false;
        }
    };
};

struct FwLn { 
    template<typename T> auto operator()(T x) -> decltype(std::log(x)) { return std::log(x); }; 
    template<typename T> bool ok(T x, Status& s) { 
        if (x>0) {
            return true;
        } else {
            s.set(Status::Domain, "Domain error.");
            return false;
        }
    }; 
};

struct FwLog10 { 
    template<typename T> auto operator()(T x) -> decltype(std::log10(x)) { return std::log10(x); }; 
    template<typename T> bool ok(T x, Status& s) { 
        if (x>0) {
            return true;
        } else {
            s.set(Status::Domain, "Domain error.");
            return false;
        }
    }; 
};

struct FwExp { 
    template<typename T> auto operator()(T x) -> decltype(std::exp(x)) { return std::exp(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true; }; 
};

struct FwSqrt { 
    template<typename T> auto operator()(T x) -> decltype(std::sqrt(x)) { return std::sqrt(x); }; 
    template<typename T> bool ok(T x, Status& s) { 
        if (x>=0) {
            return true;
        } else {
            s.set(Status::Domain, "Domain error.");
            return false;
        }
    }; 
};

struct FwAbs { 
    template<typename T> auto operator()(T x) -> decltype(std::fabs(x)) { return std::fabs(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true; }; 
};

struct FwTrunc { 
    template<typename T> auto operator()(T x) -> decltype(std::trunc(x)) { return std::trunc(x);  }; 
    template<typename T> bool ok(T x, Status& s) { return true; }; 
};

struct FwFloor { 
    template<typename T> auto operator()(T x) -> decltype(std::floor(x)) { return std::floor(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true; }; 
};

struct FwCeil { 
    template<typename T> auto operator()(T x) -> decltype(std::ceil(x)) { return std::ceil(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true; }; 
};

struct FwRound { 
    template<typename T> auto operator()(T x) -> decltype(std::round(x)) { return std::round(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true; }; 
};

struct FwSgn { 
    Real operator()(Real x) { return (x>=0) ? 1.0 : -1.0; }; 
    Int operator()(Int x) { return (x>=0) ? 1 : -1; }; 
    bool ok(Real x, Status& s) { return true; }; 
    bool ok(Int x, Status& s) { return true; }; 
};

struct FwIsInf { 
    template<typename T> Int operator()(T x) { return isInf(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true; }; 
};

struct FwIsNan { 
    template<typename T> Int operator()(T x) { return isNaN(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true; }; 
};

struct FwIsFinite { 
    template<typename T> Int operator()(T x) { return isFinite(x); }; 
    template<typename T> bool ok(T x, Status& s) { return true; }; 
};

struct FwMin { 
    Real operator()(Real x1, Real x2) { return (x1<=x2) ? x1 : x2; }; 
    Real operator()(Real x1, Int x2) { return (x1<=x2) ? x1 : x2; }; 
    Real operator()(Int x1, Real x2) { return (x1<=x2) ? x1 : x2; }; 
    Int operator()(Int x1, Int x2) { return (x1<=x2) ? x1 : x2; }; 
    bool ok(Real x1, Real x2, Status& s) { return true; }; 
};

struct FwMax { 
    Real operator()(Real x1, Real x2) { return (x1>=x2) ? x1 : x2; }; 
    Real operator()(Real x1, Int x2) { return (x1>=x2) ? x1 : x2; }; 
    Real operator()(Int x1, Real x2) { return (x1>=x2) ? x1 : x2; }; 
    Int operator()(Int x1, Int x2) { return (x1>=x2) ? x1 : x2; }; 
    bool ok(Real x1, Real x2, Status& s) { return true; }; 
};

struct FwHypot { 
    Real operator()(Real x1, Real x2) { return std::sqrt(x1*x1+x2*x2); }; 
    bool ok(Real x1, Real x2, Status& s) { return true; }; 
};

struct FwAtan2 { 
    Real operator()(Real x1, Real x2) { return std::atan(x1/x2); }; 
    bool ok(Real x1, Real x2, Status& s) { return true; }; 
};

struct FwSign { 
    Real operator()(Real x1, Real x2) { return x2>=0 ? std::fabs(x1) : -std::fabs(x1); }; 
    bool ok(Real x1, Real x2, Status& s) { return true; }; 
};

struct FwFmod { 
    Real operator()(Real x1, Real x2) { return std::fmod(x1, x2); }; 
    bool ok(Real x1, Real x2, Status& s) { 
        if (x2!=0) {
            return true;
        } else {
            s.set(Status::DivZero, "Division by zero.");
            return false;
        } 
    }; 
};

struct FwAll {
    template<typename T> Int operator()(T&& x) { 
        using Tin = typename std::remove_reference<T>::type;
        if constexpr(Value::IsVectorType<Tin>::value) {
            for(auto it=x.begin(); it!=x.end(); ++it) {
                if (!makeBool(*it)) {
                    return 0;
                }
            }
            return 1;
        } else {
            return makeBool(x);
        }
    }; 
    template<typename T> bool ok(T&& x, Status& s) { return true; }; 
};

struct FwAny {
    template<typename T> Int operator()(T&& x) { 
        using Tin = typename std::remove_reference<T>::type;
        if constexpr(Value::IsVectorType<Tin>::value) {
            for(auto it=x.begin(); it!=x.end(); ++it) {
                if (makeBool(*it)) {
                    return 1;
                }
            }
            return 0;
        } else {
            return makeBool(x);
        }
    }; 
    template<typename T> bool ok(T&& x, Status& s) { return true; }; 
};

struct FwSum {
    template<typename T> auto operator()(T&& x) -> Value::ScalarType<typename std::remove_reference<T>::type> { 
        using Tin = typename std::remove_reference<T>::type;
        if constexpr(Value::IsVectorType<Tin>::value) {
            Value::ScalarType<Tin> sum = 0;
            for(auto it=x.begin(); it!=x.end(); ++it) {
                sum += *it;
            }
            return sum;
        } else {
            return x;
        }
    }; 
    template<typename T> bool ok(T&& x, Status& s) { return true; }; 
};

struct FwProd {
    template<typename T> auto operator()(T&& x) -> Value::ScalarType<typename std::remove_reference<T>::type> { 
        using Tin = typename std::remove_reference<T>::type;
        if constexpr(Value::IsVectorType<Tin>::value) {
            Value::ScalarType<Tin> prod = 1;
            for(auto it=x.begin(); it!=x.end(); ++it) {
                prod *= *it;
            }
            return prod;
        } else {
            return x;
        }
    }; 
    template<typename T> bool ok(T&& x, Status& s) { 
        using Tin = typename std::remove_reference<T>::type;
        if constexpr(Value::IsVectorType<Tin>::value) {
            if (x.size()==0) {
                s.set(Status::Empty, "Cannot compute product of empty vector's components.");
                return false;
            }
        }
        return true; 
    }; 
};

struct FwMinAggregate {
    template<typename T> auto operator()(T&& x) -> Value::ScalarType<typename std::remove_reference<T>::type> { 
        using Tin = typename std::remove_reference<T>::type;
        if constexpr(Value::IsVectorType<Tin>::value) {
            Value::ScalarType<Tin> res = x[0];
            for(auto it=x.begin(); it!=x.end(); ++it) {
                if (*it < res) {
                    res = *it;
                }
            }
            return res;
        } else {
            return x;
        }
    }; 
    template<typename T> bool ok(T x, Status& s) { 
        using Tin = typename std::remove_reference<T>::type;
        if constexpr(Value::IsVectorType<Tin>::value) {
            Value::ScalarType<Tin> res = x[0];
            if (x.size()==0) {
                s.set(Status::Empty, "Cannot compute minimum of empty vector's components.");
                return false;
            }
        }
        return true; 
    }; 
};

struct FwMaxAggregate {
    template<typename T> auto operator()(T&& x) -> Value::ScalarType<typename std::remove_reference<T>::type> { 
        using Tin = typename std::remove_reference<T>::type;
        if constexpr(Value::IsVectorType<Tin>::value) {
            Value::ScalarType<Tin> res = 1;
            for(auto it=x.begin(); it!=x.end(); ++it) {
                if (*it > res) {
                    res = *it;
                }
            }
            return res;
        } else {
            return x;
        }
    }; 
    template<typename T> bool ok(T x, Status& s) { 
        using Tin = typename std::remove_reference<T>::type;
        if constexpr(Value::IsVectorType<Tin>::value) {
            if (x.size()==0) {
                s.set(Status::Empty, "Cannot compute maximum of empty vector's components.");
                return false;
            }
        }
        return true; 
    }; 
};

struct FwSelector {
    template<typename T> auto operator()(T& x, Int sel) -> typename T::reference {
        if (sel>=0) {
            return x[sel];
        } else {
            return x[x.size()+sel];
        }
    };
    template<typename T> bool ok(T&& x, Int sel, Status& s) {
        if (sel>=0 && sel>=x.size()) {
            s.set(Status::Range, "Selector out of range.");
            return false;
        }
        if (sel<0 && -sel>x.size()) {
            s.set(Status::Range, "Selector out of range.");
            return false;
        }
        return true;
    };
};

}

#endif
