#ifndef __ANHB_DEFINED
#define __ANHB_DEFINED

#include "parameterized.h"
#include "status.h"
#include "circuit.h"
#include "an.h"
#include "klumatrix.h"
#include "output.h"
#include "outrawfile.h"
#include "flags.h"
#include "corehb.h"
#include "common.h"


namespace NAMESPACE {

class HB : public Analysis {
public:
    typedef HBParameters Parameters;
    
    HB(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);
    
    HB           (const HB&)  = delete;
    HB           (      HB&&) = delete;
    HB& operator=(const HB&)  = delete;
    HB& operator=(      HB&&) = delete;

    virtual ~HB();
    
    virtual void dump(std::ostream& os) const;

    virtual Parameterized& parameters() { return params; }; 
    virtual const Parameterized& parameters() const { return params; }; 

    // Factory function for operating point analysis
    static Analysis* create(PTAnalysis& ptAnalysis, Circuit& circuit, Status& s=Status::ignore);

protected:
    virtual bool addCommonOutputDescriptor(const OutputDescriptor& desc);
    virtual bool addCoreOutputDescriptors(Status& s=Status::ignore);
    virtual bool resolveSave(const PTSave& save, bool verify, Status& s=Status::ignore);
    virtual bool addDefaultOutputDescriptors();
    virtual void clearOutputDescriptors();
    virtual bool resolveOutputDescriptors(bool strict, Status& s=Status::ignore);

    virtual std::tuple<bool, bool> preMapping(Status& s=Status::ignore);
    virtual bool populateStructures(Status& s=Status::ignore);

    virtual bool rebuildCores(Status& s=Status::ignore); 
    virtual bool initializeOutputs(Status& s=Status::ignore);
    virtual AnalysisCore& analysisCore() { return core; };
    virtual CoreCoroutine coreCoroutine(bool continuePrevious) {
        return std::move(core.coroutine(continuePrevious));
    };
    virtual bool formatCoreError(Status& s=Status::ignore) {
        return core.formatError(s);
    };
    virtual bool finalizeOutputs(Status& s=Status::ignore);
    virtual bool deleteOutputs(Status& s=Status::ignore);
    
    virtual size_t analysisStateStorageSize() const;
    virtual void resizeAnalysisStateStorage(size_t n);
    virtual bool storeState(size_t ndx);
    virtual bool restoreState(size_t ndx);
    virtual void makeStateIncoherent(size_t ndx);

    
private:
    IStruct<HBParameters> params;
    HBCore core;

    KluBlockSparseRealMatrix jac; // Jacobian
    VectorRepository<double> solution; // Solution history
};

}

#endif
