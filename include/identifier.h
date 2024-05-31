#ifndef __IDTABLE_DEFINED
#define __IDTABLE_DEFINED

#include <stdexcept>
#include <deque>
#include <unordered_map>
#include <string>
#include <limits>
#include <string_view>
#include <cstring>
#include <iostream>
#include "pool.h"
#include "common.h"


namespace NAMESPACE {

// Hash functor for c string
struct CstrHash {
    std::size_t operator()(const char* s) const {
        return std::hash<std::string_view>()(std::string_view(s, std::strlen(s)));
    };
};

// Comparison functors for c strings
struct CstrEqual {
    bool operator()(const char* lhs, const char* rhs) const {
        return std::strcmp(lhs, rhs)==0;
    };
};

struct CstrLess {
    bool operator()(const char* lhs, const char* rhs) const {
        return std::strcmp(lhs, rhs)<0;
    };
};

// Identifiers, common to all objects
// Everything is inline to make it fast
class Id {
public:
    // Construct id from integer
    explicit inline Id(IdentifierIndex n) : id_(n) {};

    // Default constructor constructs none identifier
    inline Id() : Id(bad) {};

    // Construct an id from string, copy string
    Id(const std::string& s);

    // Construct an id from string, move string
    Id(std::string&& s);

    // Construct an identifier from c string
    Id(const char* s);

    // Get the raw id
    inline IdentifierIndex id() const { return id_; };

    // Special 'none' id
    static const Id none;

    // Get c_str of stored string
    // inline const char* c_str() const { return identifiers[id_]->c_str(); };
    inline const char* c_str() const { 
        if (id_<nextStatic) {
            // Ordinary Id
            return indexToName()[id_]; 
        } else {
            // Static Id 
            return staticIndexToName()[maxValidIdentifierIndex-id_]; 
        }
    };

    // Convert id to string
    inline operator std::string() const { return std::string(c_str()); };

    // Convert id to bool
    inline operator bool() const { return id_ != bad; };
    
    // Compare two identifiers
    inline bool operator==(const Id& other) const { return id_ == other.id_; };
    inline bool operator!=(const Id& other) const { return id_ != other.id_; };

    // Create static identifier
    // Static identifiers should be created as global variables with static linkage. 
    // This way the static identifiers are created before any ordinary identifier. 
    // Static identifiers should never be registered in dynamically loaded modules
    // because their place may already be taken by an existing ordinary identifier. 
    static Id createStatic(const char* s);
    
    // Write to stream
    friend std::ostream& operator<<(std::ostream& os, const Id& obj);

    static void dump(std::ostream& os, bool details=false) { stringPool().dump(os, details); };

private:
    // 2**32-1 possible identifiers, MAX reserved for Id::none
    static const IdentifierIndex maxValidIdentifierIndex = std::numeric_limits<IdentifierIndex>::max();
    static const IdentifierIndex bad = 0;

    void constructorHelper(const char *s, size_t n, bool makeStatic=false);

    // Actual id
    IdentifierIndex id_;

    // Pool of c strings
    static CStringPool& stringPool();
    static CStringPool& staticStringPool();

    // Map from const char* to index
    static std::unordered_map<const char*,IdentifierIndex,CstrHash,CstrEqual>& nameToIndex();
    
    // Map from index to const char*
    static std::deque<const char*>& indexToName();
    static std::deque<const char*>& staticIndexToName();
    
    static IdentifierIndex nextOrdinary;
    static IdentifierIndex nextStatic;

#ifdef SIMDEBUG
    const char* str; 
#endif
};

}

// Do this in std namespace because there are many maps using Id as key and 
// we don't want to explicitnly specify the hash function each time. 
namespace std {
    // Hash functor for Id (specialization of std::hash)
    template<> struct hash<NAMESPACE::Id> {
        std::size_t operator()(const NAMESPACE::Id& k) const {
            return hash<NAMESPACE::IdentifierIndex>()(k.id());
        };
    };
}


#endif
