#include "libplatform.h"
#include "simulator.h"
#include "parser.h"
#include "openvafcomp.h"
#include "circuit.h"
#include "processutils.h"

using namespace sim;

// Circuit building and analysis demo (with constant options and saves)
int main() {
    // Path to staged models (osdi files)
    std::string modulePath = "../../lib/vacask/mod";
    // Path to staged Python libraries
    std::string pythonLibraryPath = "../../lib/vacask/python";
    // Python binary name
    std::string pythonBinary = findPythonExecutable();

    // Status variable
    Status s;

    // Simulator setup, no paths set
    Simulator::setup();

    // Prepend directories to the list of osdi file paths
    Simulator::prependModulePath({modulePath});

    // Parser tables
    ParserTables tab("RC transient");

    // Parser, needs tables to store stats when parsing expressions and parameters
    Parser p(tab);

    // Title, loads, default ground (0)
    tab
        .add(PTLoad("resistor.osdi"))
        .add(PTLoad("capacitor.osdi"))
        .defaultGround()
        .setDefaultSubDef(
            // Toplevel subcircuit definition, no name, no terminals
            PTSubcircuitDefinition()
            .add(PTParameters(PVv( PV{"c0", 1e-6}, PV{"v0", 5} ), {}))
            // Resistor model res, capacitor model cap, voltage source model vsrc
            .add(PTModel("res", "resistor"))
            .add(PTModel("cap", "capacitor"))
            .add(PTModel("vsrc", "vsource"))
            // Element r1, master is res, terminals 1 and 2, constant parameter r=1000, no parameter expressions
            .add(PTInstance("r1", "res", {"1", "2"})
                // Add constant parameter
                .add(PV{"r", 1000})
            )
            // Element c1, master cap, terminals 1 and 0, constant parameter c=1e-6, no parameter expressions
            .add(PTInstance("c1", "cap", {"2", "0"})
                // Add parameter defined with an expression
                .add(PE{"c", p.parseExpression("2*c0")})
            )
            // Pulse source
            .add(PTInstance("v1", "vsrc", {"1", "0"})
                // Add parsed parameters
                .add(p.parseParameters("type=\"pulse\" val0=0 val1=v0 delay=1m rise=1u fall=1u width=4m"))
            )
        )
        // Embedded Python postprocessing script
        .add(PTEmbed("runme.py", R"script(
from rawfile import rawread
import numpy as np
import matplotlib.pyplot as plt 
import sys

tran1 = rawread('tran1.raw').get()
print("Vectors:", tran1.names)
fig1, ax1 = plt.subplots(1, 1, figsize=(6,4), dpi=100, constrained_layout=True)
fig1.axes[0].plot(tran1["time"], tran1["2"], "r", marker=".")
fig1.axes[0].plot(tran1["time"], tran1["1"], "b")
fig1.axes[0].plot(tran1["time"], tran1["r1.i"]*1000, "--")

plt.show()
)script"));

    // Verify tables
    if (!tab.verify(s)) {
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    // Dump tables for debugging
    tab.dump(0, Simulator::out());

    // Store embedded files (we don't have any, but this is how you do it)
    if (!tab.writeEmbedded(1, s)) {
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    // Create circuit, create OpenVAF compiler with no options
    OpenvafCompiler comp;
    // Circuit object
    Circuit cir(tab, &comp, s);
    if (!cir.isValid()) {
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    // Change reltol
    cir.setOption("reltol", 1e-4);
    
    // You could also do
    // IStruct<SimulatorOptions> opt;
    // opt.core().reltol = 1e-4;
    // cir.setOptions(opt);

    // Do not run this in a loop, stuff shoud be parsed only once 
    // because all strings are stored in tab. 
    // auto parsedOpt = p.parseParameters("reltol=1e-4");
    // cir.setOptions(parsedOpt);
    
    // Elaborate it, just the default toplevel subcircuit 
    // (empty list of subcircuit definition names). 
    // Use __topdef__ and __topinst__ as prefixes for toplevel 
    // subcircuit model and toplevel subcircuit instance names. 
    // Do not collect device requests (i.e. abort/finish/stop). 
    if (!cir.elaborate({}, "__topdef__", "__topinst__", nullptr, s)) {
        Simulator::err() << "Elaboration failed.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    cir.dumpHierarchy(0, Simulator::out());

    // Analysis description
    auto tranDesc = PTAnalysis("tran1", "tran");
    tranDesc
        .add(PV{"step", 1e-6})
        .add(PV("stop", 10e-3));

    // Save directives
    
    // Analysis object, no options map (options that are given as parameterized expressions)
    auto tran = Analysis::create(tranDesc, cir, s);
    tran->add(PTSave("default"));
    tran->add(PTSave("p", "r1", "i"));
    if (!tran) {
        Simulator::err() << "Failed to create analysis.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    // Run analysis
    auto [ok, canResume] = tran->run(s);
    if (!ok) {
        Simulator::err() << "Analysis failed.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }
    Simulator::out() << "Analysis OK. Can resume: " << (canResume ? "true" : "false") << "\n";

    // Cleanup
    delete tran;

    // Run postprocessing
    runProcess(pythonBinary, {"runme.py"}, &pythonLibraryPath, false, false);

    return 0;
}
