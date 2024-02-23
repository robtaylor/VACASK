#ifndef __ANACTF_DEFINED
#define __ANACTF_DEFINED

#include "an.h"
#include "coreop.h"
#include "coreactf.h"
#include "parameterized.h"
#include "common.h"


namespace NAMESPACE {

class AcTf : public Analysis {
public:
    typedef AcTfParameters Parameters;

    AcTf(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);
    
    AcTf           (const AcTf&)  = delete;
    AcTf           (      AcTf&&) = delete;
    AcTf& operator=(const AcTf&)  = delete;
    AcTf& operator=(      AcTf&&) = delete;

    virtual ~AcTf();

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
    IStruct<AcTfParameters> params;
    OperatingPointCore opCore;
    AcTfCore acTfCore;
    std::unordered_map<Id,size_t> sourceIndex;
    
    KluRealMatrix jac; // Resistive Jacobian
    VectorRepository<double> solution; // Solution history
    VectorRepository<double> states; // Circuit states

    KluComplexMatrix acMatrix; 
    Vector<Complex> acSolution;

    std::vector<Instance*> sources;
    Vector<Complex> tf;
    Vector<Complex> yin;
    Vector<Complex> zin;
};

}

#endif
