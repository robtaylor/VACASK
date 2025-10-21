#ifndef __NATURES_DEFINED
#define __NATURES_DEFINED

#include <string>
#include <map>
#include <unordered_set>
#include "identifier.h"
#include "common.h"


namespace NAMESPACE {

// Used by devices for returning node/residual tolerances. 
// When device is asked for tolerances, it returns tolerance and nature id. 
// At tolerance setup id and tolerance are stored in the corresponding vector. 
// If a tolerance is updated, the node's nature changes to the latest nature. 
// After tolerances are collected ids are converted to numbers via registry. 
// Nonlinear solvers work with nature numbers. 
// Nature numbers are stored in CommonData vectors. 
// They are used for setting up global maxima (point and historic). 
using NatureId = Id;

typedef struct NatureTolerance {
    NatureId id;
    double abstol;
    NatureTolerance(NatureId id, double abstol) : id(id), abstol(abstol) {};
} NatureTolerance;


// Natures are primarily identified by name. 
// For ordinary natures the identifier is based on the nature <name>.
// For disciplines the identifier is based on <discipline>.potential or <discipline>.flow.
// Builtin (spice) natures have identifiers based on .<name> (.voltage, .current, ...).
// Stores all files where a nature is found.
// Issues nature identifiers.
// There is one such global registry (static members). 
// Each file registers natures and obtains their ids when its dynamic library is first opened. 
class NatureRegistry {
public:
    NatureRegistry           (const NatureRegistry&)  = delete;
    NatureRegistry           (      NatureRegistry&&) = delete;
    NatureRegistry& operator=(const NatureRegistry&)  = delete;
    NatureRegistry& operator=(      NatureRegistry&&) = delete;

    // Obtain nature Id and add file to the nature's list of files
    static NatureId natureId(std::string file, std::string name);

    // NatureId values for builtin natures (no nature and SPICE natures)
    static const NatureId noNature;
    static const NatureId spiceVoltage;
    static const NatureId spiceCurrent;
    static const NatureId spiceFlux;
    static const NatureId spiceCharge;

private:
    static std::unordered_map<NatureId,std::unordered_set<std::string>> natureFilesMap;

    static NatureId natureIdInternal(std::string file, std::string name, bool staticId=false);
};


// Subset of natures used by a circuit. 
// Maps nature identifier to index and vice versa. 
// There is one such object per circuit. 
// Each device registers its natures (ids) when circuit's loads are processed. 
class NaturesSubset {
public:
    using NatureIndex = uint32_t;

    NaturesSubset           (const NaturesSubset&)  = delete;
    NaturesSubset           (      NaturesSubset&&) = delete;
    NaturesSubset& operator=(const NaturesSubset&)  = delete;
    NaturesSubset& operator=(      NaturesSubset&&) = delete;

    NaturesSubset();

    // Reset
    void reset();

    // Number of registered natures
    NatureIndex count() const { return indexToId.size(); };

    // Add nature to registry
    NatureIndex registerNature(NatureId id);

    // Obtain nature's index
    NatureIndex natureIndex(NatureId id) { 
        if (idToIndex.contains(id)) {
            return idToIndex[id];
        } else {
            return NaturesSubset::noNature;
        }
    };

    // Obtain nature Id from index
    // Does not perform range check. 
    NatureId natureId(NatureIndex ndx) { return indexToId[ndx]; };

    // Indices of no nature and SPICE natures
    // No nature is index 0. 
    static const NatureIndex noNature;
    static const NatureIndex spiceVoltage;
    static const NatureIndex spiceCurrent;
    static const NatureIndex spiceFlux;
    static const NatureIndex spiceCharge;

private:
    std::unordered_map<NatureId,NatureIndex> idToIndex;
    std::vector<NatureId> indexToId;
};

}

#endif
