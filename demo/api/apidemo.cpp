#include "sourceloc.h"
#include "parseroutput.h"
#include "simulator.h"
#include "circuit.h"
#include "an.h"
#include "processutils.h"


using namespace sim;

// TODO: 
//   Common linking into executable for all binaries
//   API functions without Loc
//   API functions for chaining, i.e. add()
//   Move parser in lib
//   String parsing into 
//     RPN
//     PTParameters
int main() {
    // Path to staged models (osdi files)
    std::string modulePath = "../../lib/vacask/mod";
    // Path to staged Python libraries
    std::string pythonLibraryPath = "../../lib/vacask/python";
    // Python binary name
    std::string pythonBinary = "python3";

    // Status variable
    Status s;

    // Simulator setup, no paths set
    Simulator::setup();

    // Prepend directories to the list of osdi file paths
    Simulator::prependModulePath({modulePath});

    // Parser tables
    ParserTables tab("RC transient");

    // No location, we are not parsing an input file or a string
    auto nl = Loc::bad;

    // Title, loads, default ground (0)
    tab
        .addLoad(PTLoad(nl, "resistor.osdi"))
        .addLoad(PTLoad(nl, "capacitor.osdi"))
        .defaultGround();

    // Toplevel subcircuit definition named __topdef__, no terminals
    auto toplevel = PTSubcircuitDefinition(nl, "__topdef__", {});
    // Get root block, add models and instances
    toplevel.root()
        // Resistor model res, capacitor model cap, voltage source model vsrc
        .add(PTModel(nl, "res", "resistor"))
        .add(PTModel(nl, "cap", "capacitor"))
        .add(PTModel(nl, "vsrc", "vsource"))
        // Element r1, master is res, terminals 1 and 2, constant parameter r=1000, no parameter expressions
        .add(PTInstance(nl, "r1", "res", {"1", "2"}, PTParameters(PVv( PV{nl, "r", 1000} ), {})))
        // Element c1, master cap, terminals 1 and 0, constant parameter c=1e-6, no parameter expressions
        .add(PTInstance(nl, "c1", "cap", {"2", "0"}, PTParameters(PVv( PV{nl, "c", 1e-6} ), {})))
        // Pulse source
        .add(PTInstance(nl, "v1", "vsrc", {"1", "0"}, PTParameters(PVv( 
            PV{nl, "type", "pulse"}, 
            PV{nl, "val0", 0}, 
            PV{nl, "val1", 5}, 
            PV{nl, "delay", 1e-3}, 
            PV{nl, "rise", 1e-6}, 
            PV{nl, "fall", 1e-6}, 
            PV{nl, "width", 4e-3} 
        ), {})));

    // Move toplevel subcircuit definition into ParserTables as default subdef
    tab.addDefaultSubDef(std::move(toplevel));

    // Dump tables for debugging
    tab.dump(0, Simulator::out());

    // Create circuit, no SourceCompiler
    Circuit cir(tab, nullptr, s);
    if (!cir.isValid()) {
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    // Elaborate it, just the default toplevel subcircuit. 
    // Use __topdef__ and __topinst__ as prefixes for toplevel subcircuit definition and instance names. 
    // No options, do not collect device requests. 
    if (!cir.elaborate({}, "__topdef__", "__topinst__", nullptr, nullptr, s)) {
        Simulator::err() << "Elaboration failed.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    cir.dumpHierarchy(0, Simulator::out());

    // Analysis description
    auto tranDesc = PTAnalysis(nl, "tran1", "tran");
    tranDesc.add(PTParameters(PVv(
        PV{nl, "step", 1e-6}, 
        PV(nl, "stop", 10e-3)
    ), {}));

    // Analysis object, no saves, no options map (options that are given as parameterized expressions)
    auto tran = Analysis::create(tranDesc, nullptr, nullptr, cir, s);
    if (!tran) {
        Simulator::err() << "Failed to create analysis.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }
    auto [ok, canResume] = tran->run(s);
    if (!ok) {
        Simulator::err() << "Analysis failed.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }
    Simulator::out() << "Analysis OK. Can resume: " << (canResume ? "true" : "false") << "\n";

    // Cleanup
    delete tran;

    // Python postprocessing script
    std::ofstream file("runme.py");
    file << R"script(
from rawfile import rawread
import numpy as np
import matplotlib.pyplot as plt 
import sys

tran1 = rawread('tran1.raw').get()
fig1, ax1 = plt.subplots(1, 1, figsize=(6,4), dpi=100, constrained_layout=True)
fig1.axes[0].plot(tran1["time"], tran1["2"], "r", marker=".")
fig1.axes[0].plot(tran1["time"], tran1["1"], "b")

plt.show()
)script";
    file.close();
    runProcess(pythonBinary, {"runme.py"}, &pythonLibraryPath, false, false);

    return 0;
}
