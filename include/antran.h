#ifndef __ANTRAN_DEFINED
#define __ANTRAN_DEFINED

#include "an.h"
#include "coreop.h"
#include "coretran.h"
#include "parameterized.h"
#include "common.h"


namespace NAMESPACE {

class Tran : public Analysis {
public:
    typedef TranParameters Parameters;

    Tran(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);
    
    Tran           (const Tran&)  = delete;
    Tran           (      Tran&&) = delete;
    Tran& operator=(const Tran&)  = delete;
    Tran& operator=(      Tran&&) = delete;

    virtual ~Tran();

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
    virtual AnalysisCore& analysisCore() { return tranCore; };
    virtual CoreCoroutine coreCoroutine(bool continuePrevious) {
        return std::move(tranCore.coroutine(continuePrevious));
    };
    virtual bool formatCoreError(Status& s=Status::ignore) {
        return tranCore.formatError(s);
    };
    virtual bool finalizeOutputs(Status& s=Status::ignore);
    virtual bool deleteOutputs(Status& s=Status::ignore);
    
    virtual size_t analysisStateStorageSize() const;
    virtual size_t allocateAnalysisStateStorage(size_t n);
    virtual void deallocateAnalysisStateStorage(size_t n=0);
    virtual bool storeState(size_t ndx, bool storeDetails=true);
    virtual bool restoreState(size_t ndx);
    virtual void makeStateIncoherent(size_t ndx);

    
private:
    IStruct<TranParameters> params;
    OperatingPointCore opCore;
    TranCore tranCore;
    
    KluRealMatrix jac; // Jacobian
    VectorRepository<double> solution; // Solution history
    VectorRepository<double> states; // Circuit states
};

}

#endif
