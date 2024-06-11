#ifndef __ACCT_DEFINED
#define __ACCT_DEFINED

#include "common.h"
#include <chrono>

namespace NAMESPACE {

typedef struct AcctData {
    int parse {0};
    double tparse {0};

    int elab {0};
    double telab {0};

    int chgelab {0};
    double tchgelab {0};

    int load {0};
    int resload {0};
    int sysload {0};
    int factor {0};
    int refactor {0};
    int solve {0};
    int cxfactor {0};
    int cxrefactor {0};
    int cxsolve {0};
    int nrcall {0};
    int nriter {0};
    int accepted {0};
    int rejected {0};
    int points {0};
    double tload {0};
    double tresload {0};
    double tsysload {0};
    double tfactor {0};
    double trefactor {0}; 
    double tsolve {0};
    double tcxfactor {0};
    double tcxrefactor {0}; 
    double tcxsolve {0};
    double tnr {0};
    double t {0};
} AcctData;

typedef struct Accounting {
    AcctData acctOld;
    AcctData acctNew;
    
    /*void advanceParse() { acctOld.parse = acctNew.parse; };
    void advanceAnalysis() { acctOld.analysis = acctNew.analysis; };
    void advanceSweep() { acctOld.sweep = acctNew.sweep; };
    void advancePoint() { acctOld.point = acctNew.point; };*/

    // High resolution wall clock
    typedef decltype(std::chrono::high_resolution_clock::now()) Timepoint; 
    static Timepoint wclk() { return std::chrono::high_resolution_clock::now(); }; 
    static double wclkDelta(Timepoint t0) { 
        auto interval = std::chrono::high_resolution_clock::now() - t0;
        return std::chrono::duration_cast<std::chrono::duration<double>>(interval).count();
    }; 

    void dumpTotal(int indent, std::ostream& os) const;
    void dumpChange(int indent, std::ostream& os) const;
} Accounting;

}

#endif
