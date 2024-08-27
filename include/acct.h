#ifndef __ACCT_DEFINED
#define __ACCT_DEFINED

#include "common.h"
#include <chrono>

namespace NAMESPACE {

typedef struct AcctData {
    size_t parse {0};
    double tparse {0};

    size_t elab {0};
    double telab {0};

    size_t chgelab {0};
    double tchgelab {0};

    size_t conv {0};
    size_t evalload {0};
    size_t resload {0};
    size_t sysload {0};
    size_t factor {0};
    size_t refactor {0};
    size_t solve {0};
    size_t cxfactor {0};
    size_t cxrefactor {0};
    size_t cxsolve {0};
    size_t nrcall {0};
    size_t nriter {0};
    size_t accepted {0};
    size_t rejected {0};
    size_t points {0};
    size_t bpiicount {0};
    size_t bpiibypass {0};
    size_t bpiibpfailed {0};
    size_t bpiiconvcheck {0};
    size_t bpiiconverged {0};
    double tconv {0};
    double tevalload {0};
    double teval {0};
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
    AcctData acctPrevParse;
    AcctData acctPrevAnalysis;
    AcctData acctPrevSweepPoint;
    AcctData acctPrevPoint;
    AcctData acctNew;
    
    // High resolution wall clock
    typedef decltype(std::chrono::high_resolution_clock::now()) Timepoint; 
    static Timepoint wclk() { return std::chrono::high_resolution_clock::now(); }; 
    static double timeDelta(Timepoint& begin, Timepoint& end) { 
        return std::chrono::duration_cast<std::chrono::duration<double>>(end-begin).count(); 
    };
    static double wclkDelta(Timepoint& t0, bool update=false) { 
        auto now = std::chrono::high_resolution_clock::now();
        if (update) {
            t0 = now;
        }
        return timeDelta(t0, now);
    }; 

    void dumpTotal(int indent, std::ostream& os) const;
    void dumpChange(int indent, std::ostream& os) const;
} Accounting;

}

#endif
