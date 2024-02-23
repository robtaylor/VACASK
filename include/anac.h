#ifndef __ANAC_DEFINED
#define __ANAC_DEFINED

#include "an.h"
#include "coreop.h"
#include "coreac.h"
#include "parameterized.h"
#include "common.h"


namespace NAMESPACE {

class Ac : public Analysis {
public:
    typedef AcParameters Parameters;

    Ac(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);
    
    Ac           (const Ac&)  = delete;
    Ac           (      Ac&&) = delete;
    Ac& operator=(const Ac&)  = delete;
    Ac& operator=(      Ac&&) = delete;

    virtual ~Ac();

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
    IStruct<AcParameters> params;
    OperatingPointCore opCore;
    AcCore acCore;
    
    KluRealMatrix jac; // Resistive Jacobian
    VectorRepository<double> solution; // Solution history
    VectorRepository<double> states; // Circuit states

    KluComplexMatrix acMatrix; 
    Vector<Complex> acSolution;
};

}

#endif
