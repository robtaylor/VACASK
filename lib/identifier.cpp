#include "identifier.h"
#include <iostream>
#include "common.h"


namespace NAMESPACE {

// Pool 

CStringPool& Id::stringPool() {
    static CStringPool internal(stringPoolBlockSize, stringPoolGrowthFactor, stringPoolRetries);
    return internal;
}
CStringPool& Id::staticStringPool() {
    static CStringPool internal(stringPoolBlockSize, stringPoolGrowthFactor, stringPoolRetries);
    return internal;
}

std::unordered_map<const char*,IdentifierIndex,CstrHash,CstrEqual>& Id::nameToIndex() {
    static std::unordered_map<const char*,IdentifierIndex,CstrHash,CstrEqual> internal;
    return internal;
}

std::deque<const char*>& Id::indexToName() {
    static const char* badIdCstring ="";
    static std::deque<const char*> internal({badIdCstring});
    return internal;
}

std::deque<const char*>& Id::staticIndexToName() {
    static std::deque<const char*> internal({});
    return internal;
}

IdentifierIndex Id::nextOrdinary = 1;
IdentifierIndex Id::nextStatic = Id::maxValidIdentifierIndex; 

void Id::constructorHelper(const char *s, size_t n, bool makeStatic) {
    // Look it up
    auto it = nameToIndex().find(s); 
    if (std::string(s)=="instance") {
        auto a=1;
    }
    if (it!=nameToIndex().end()) {
        // Exists
        id_ = it->second;
        if (makeStatic && id_<nextOrdinary) {
            // Ordinary identifier with same name already exists
            throw std::logic_error("Static identifier conflicts an ordinary identifier ("+std::string(s)+").");
        }
   } else {
        // Get new index
        IdentifierIndex i;
        const char* cstr;
        if (makeStatic) {
            // Static identifier
            i = nextStatic;
            if (i<=nextOrdinary) {
                throw std::length_error("Too many identifiers.");
            }
            nextStatic--;

            // Put into pool
            cstr = staticStringPool().allocate(s, n);

            // Add to staticIndexToName
            staticIndexToName().push_back(cstr);
        } else {
            // Ordinary identifier
            i = nextOrdinary;
            if (i>=nextStatic) {
                throw std::length_error("Too many identifiers.");
            }
            nextOrdinary++;

            // Put into pool
            cstr = stringPool().allocate(s, n);

            // Add to indexToName
            indexToName().push_back(cstr);
        }
        
        // Insertinto nameToIndex works always because find() found nothing
        bool inserted;
        std::tie(it, inserted) = nameToIndex().insert({cstr, i});
        
        id_ = i;
    }
#ifdef SIMDEBUG
    str = it->first;
#endif
    // std::cout << "copy " << id_ << " " << it->first << "\n"";
}

Id::Id(const std::string& s) {
    constructorHelper(s.c_str(), s.size());
}

Id::Id(std::string&& s) {
    constructorHelper(s.c_str(), s.size());
}

Id::Id(const char* s) {
    constructorHelper(s, std::strlen(s));
}

Id Id::createStatic(const char* s) {
    auto id = Id();
    id.constructorHelper(s, std::strlen(s), true);
    return id;
}

const Id Id::none = Id();

std::ostream& operator<<(std::ostream& os, const Id& id) {
    os << std::string(id);
    return os;
}

}
