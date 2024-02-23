#ifndef __ANDCTF_DEFINED
#define __ANDCTF_DEFINED

#include "an.h"
#include "coreop.h"
#include "coredctf.h"
#include "parameterized.h"
#include "common.h"


namespace NAMESPACE {

class DcTf : public Analysis {
public:
    typedef DcTfParameters Parameters;

    DcTf(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);
    
    DcTf           (const DcTf&)  = delete;
    DcTf           (      DcTf&&) = delete;
    DcTf& operator=(const DcTf&)  = delete;
    DcTf& operator=(      DcTf&&) = delete;

    virtual ~DcTf();

    virtual void dump(std::ostream& os) const;

    virtual Parameterized& parameters() { return params; }; 
    virtual const Parameterized& parameters() const { return params; }; 

    // Factory function for operating point analysis
    static Analysis* create(PTAnalysis& ptAnalysis, Circuit& circuit, Status& s=Status::ignore);

protected:
    virtual bool addCommonOutputDescriptor(const OutputDescriptor& desc);
    virtual bool addCoreOutputDescriptors(Status& s=Status::ignore);
    virtual bool resolveSave(const PTSave& save, bool verify, Status& s=Status::ignore);
    virtual bool addDefaultOutputDescriptors(Status& s=Status::ignore);
    virtual bool resolveOutputDescriptors(bool strict, Status& s=Status::ignore);

    virtual std::tuple<bool, bool> preMapping(Status& s=Status::ignore);
    virtual bool populateStructures(Status& s=Status::ignore);

    virtual bool rebuildCores(Status& s=Status::ignore); 
    virtual bool initializeOutputs(Status& s=Status::ignore);
    virtual bool runCores(bool continuePrevious, Status& s=Status::ignore);
    virtual bool finalizeOutputs(Status& s=Status::ignore);
    virtual bool deleteOutputs(Status& s=Status::ignore);
    
    virtual size_t analysisStateStorageSize() const;
    virtual void resizeAnalysisStateStorage(size_t n);
    virtual bool storeState(size_t ndx);
    virtual bool restoreState(size_t ndx);
    virtual void makeStateIncoherent(size_t ndx);

    virtual void clearOutputDescriptors();
    
private:
    IStruct<DcTfParameters> params;
    OperatingPointCore opCore;
    DcTfCore dcTfCore;
    std::unordered_map<Id,size_t> sourceIndex;
    std::vector<Instance*> sources; 
    
    KluRealMatrix jac; // Resistive Jacobian
    VectorRepository<double> solution; // Solution history
    VectorRepository<double> states; // Circuit states

    Vector<double> incrementalSolution; // Incremental solution
    Vector<double> tf;
    Vector<double> yin;
    Vector<double> zin;

};

}

#endif
