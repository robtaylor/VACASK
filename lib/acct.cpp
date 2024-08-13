#include "acct.h"
#include "common.h"

namespace NAMESPACE {

void Accounting::dumpTotal(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "Parser invocations:              " << acctNew.parse << "\n";
    os << pfx << "Parse time:                      " << acctNew.tparse << "\n";

    os << "\n";

    os << pfx << "Elaborations:                    " << acctNew.elab << "\n";
    os << pfx << "Elaboration time:                " << acctNew.telab << "\n";

    os << "\n";

    os << pfx << "Partial elaborations:            " << acctNew.chgelab << "\n";
    os << pfx << "Partial elaboration time:        " << acctNew.tchgelab << "\n";

    os << "\n";

    os << pfx << "Eval and load calls:             " << acctNew.evalload << "\n";
    os << pfx << "Eval and load time:              " << acctNew.tevalload << "\n";
    // os << pfx << "Load time:                      " << acctNew.tload << "\n";
    // os << pfx << "Jacobian write time:            " << acctNew.teval << "\n";

    os << "\n";

    os << pfx << "Instance convergence checks:     " << acctNew.conv << "\n";
    os << pfx << "Instance convergence check time: " << acctNew.tconv << "\n";

    os << "\n";

    os << pfx << "NR solver calls:                 " << acctNew.nrcall << "\n";
    os << pfx << "NR solver iterations:            " << acctNew.nriter << "\n";
    os << pfx << "NR solver time:                  " << acctNew.tnr << "\n";

    os << "\n";

    if (acctNew.bpiicount>0) {
        os << pfx << "Bypassable instance iterations:  " << acctNew.bpiicount << "\n";
        os << pfx << "Converged:                       " << acctNew.bpiiconv << "\n";
        os << pfx << "Bypassed:                        " << acctNew.bpiibypass << "\n";
        os << pfx << "Bypasses failed:                 " << acctNew.bpiibpfailed << "\n";
        os << "\n";
    }

    os << pfx << "Real factorizations:             " << acctNew.factor << "\n";
    os << pfx << "Real refactorizations:           " << acctNew.refactor << "\n";
    os << pfx << "Real solve calls:                " << acctNew.solve << "\n";
    os << pfx << "Real factorization time:         " << acctNew.tfactor << "\n";
    os << pfx << "Real refactorization time:       " << acctNew.trefactor << "\n";
    os << pfx << "Real solve time:                 " << acctNew.tsolve << "\n";
    
    os << "\n";

    os << pfx << "Complex factorizations:          " << acctNew.cxfactor << "\n";
    os << pfx << "Complex refactorizations:        " << acctNew.cxrefactor << "\n";
    os << pfx << "Complex solve calls:             " << acctNew.cxsolve << "\n";
    os << pfx << "Complex factorization time:      " << acctNew.tcxfactor << "\n";
    os << pfx << "Complex refactorization time:    " << acctNew.tcxrefactor << "\n";
    os << pfx << "Complex solve time:              " << acctNew.tcxsolve << "\n";

    os << "\n";

    os << pfx << "Accepted timepoints:             " << acctNew.accepted << "\n";
    os << pfx << "Rejected timepoints:             " << acctNew.rejected << "\n";
}

}
