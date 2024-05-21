#ifndef __ANOP_DEFINED
#define __ANOP_DEFINED

#include "parameterized.h"
#include "status.h"
#include "circuit.h"
#include "an.h"
#include "klumatrix.h"
#include "output.h"
#include "outrawfile.h"
#include "flags.h"
#include "coreop.h"
#include "common.h"


namespace NAMESPACE {

class OperatingPoint : public Analysis {
public:
    typedef OpParameters Parameters;
    
    OperatingPoint(Id name, Circuit& circuit, PTAnalysis& ptAnalysis);
    
    OperatingPoint           (const OperatingPoint&)  = delete;
    OperatingPoint           (      OperatingPoint&&) = delete;
    OperatingPoint& operator=(const OperatingPoint&)  = delete;
    OperatingPoint& operator=(      OperatingPoint&&) = delete;

    virtual ~OperatingPoint();
    
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
    virtual bool runCores(bool continuePrevious, Status& s=Status::ignore);
    virtual bool finalizeOutputs(Status& s=Status::ignore);
    virtual bool deleteOutputs(Status& s=Status::ignore);
    
    virtual size_t analysisStateStorageSize() const;
    virtual void resizeAnalysisStateStorage(size_t n);
    virtual bool storeState(size_t ndx);
    virtual bool restoreState(size_t ndx);
    virtual void makeStateIncoherent(size_t ndx);

    
private:
    IStruct<OpParameters> params;
    OperatingPointCore core;

    KluRealMatrix jac; // Resistive Jacobian
    VectorRepository<double> solution; // Solution history
    VectorRepository<double> states; // Circuit states
};

}

#endif
