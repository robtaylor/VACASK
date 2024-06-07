#include "acct.h"
#include "common.h"

namespace NAMESPACE {

void Accounting::dumpTotal(int indent, std::ostream& os) const {
    std::string pfx = std::string(indent, ' ');
    os << pfx << "Parser invocations:       " << acctNew.parse.parse << "\n";
    os << pfx << "Parse time:               " << acctNew.parse.tparse << "\n";

    os << "\n";

    os << pfx << "Elaborations:             " << acctNew.analysis.elab << "\n";
    os << pfx << "Elaboration time:         " << acctNew.analysis.telab << "\n";

    os << "\n";

    os << pfx << "Partial elaborations:     " << acctNew.sweep.chgelab << "\n";
    os << pfx << "Partial elaboration time: " << acctNew.sweep.tchgelab << "\n";

    os << "\n";

    os << pfx << "Eval and load calls:      " << acctNew.point.load << "\n";
    os << pfx << "Eval and load time:       " << acctNew.point.tload << "\n";

    os << "\n";

    os << pfx << "NR solver calls:          " << acctNew.point.nrcall << "\n";
    os << pfx << "NR solver iterations:     " << acctNew.point.nriter << "\n";
    os << pfx << "NR solver time:           " << acctNew.point.tnr << "\n";

    os << "\n";

    os << pfx << "Accepted timepoints:      " << acctNew.point.accepted << "\n";
    os << pfx << "Rejected timepoints:      " << acctNew.point.rejected << "\n";

}

}
