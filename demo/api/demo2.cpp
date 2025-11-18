#include "libplatform.h"
#include "simulator.h"
#include "parser.h"
#include "openvafcomp.h"
#include "circuit.h"
#include "processutils.h"

using namespace sim;

// Parameterized subcircuit with conditional blocks and a DC sweep
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
    ParserTables tab("Conditional blocks");

    // Parser, needs tables to store stats when parsing expressions and parameters
    Parser p(tab);

    // Title, loads, default ground (0)
    tab
        .add(PTLoad("resistor.osdi"))
        .defaultGround()
        .setDefaultSubDef(
            // Toplevel subcircuit definition, no name, no terminals
            PTSubcircuitDefinition()
            // Resistor model res, capacitor model cap, voltage source model vsrc
            .add(PTModel("res", "resistor"))
            .add(PTModel("vsrc", "vsource"))
            // Add mysub subcircuit
            .add(PTSubcircuitDefinition("mysub", { "p", "out", "n" })
                .add(PV("r", 1000))
                .add(PV("fact", 0.5))
                .add(PV("mode", 0))
                .add(PTBlockSequence()
                    .add(p.parseExpression("mode==0"), PTBlock()
                        .add(PTInstance("r1", "res", {"p", "out"})
                            .add(PE("r", p.parseExpression("r*(1-fact)")))
                        )
                        .add(PTInstance("r2", "res", {"out", "n"})
                            .add(PE("r", p.parseExpression("r*fact")))
                        )
                    )
                    // Else block has an empty Rpn object as condition
                    .add(Rpn(), PTBlock()
                        .add(PTInstance("r1", "res", {"p", "out"})
                            .add(PE("r", p.parseExpression("r*fact")))
                        )
                        .add(PTInstance("r2", "res", {"out", "n"})
                            .add(PE("r", p.parseExpression("r*(1-fact)")))
                        )
                    )
                )
            )
            // Instance of mysub
            .add(PTInstance("x1", "mysub", {"1", "2", "0"})
                .add(PV{"r", 1000})
                .add(PV{"fact", 0.2})
            )
            // Voltage source
            .add(PTInstance("v1", "vsrc", {"1", "0"})
                .add(PV{"dc", 10})
            )
        )
        // Embedded Python postprocessing script
        .add(PTEmbed("runme.py", R"script(
from rawfile import rawread
import numpy as np
import matplotlib.pyplot as plt 
import sys

dc1 = rawread('dc1.raw').get()
print("Vectors:", dc1.names)
fig1, ax1 = plt.subplots(1, 1, figsize=(6,4), dpi=100, constrained_layout=True)
print("Mode:", dc1["mode"])
print("Output:", dc1["2"])
print("Expected: [2, 8]")
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

    // Set circuit variable var1 to 60
    cir.setVariable("var1", 60);

    // Elaborate it, just the default toplevel subcircuit 
    // (empty list of subcircuit definition names). 
    // Use __topdef__ and __topinst__ as prefixes for toplevel 
    // subcircuit model and toplevel subcircuit instance names. 
    // Do not collect device requests (i.e. abort/finish/stop). 
    // This time we use an initializer list of PTParameters objects for initial simulator options
    // The first PTParameters object holds the reltol value. 
    // The second one holds an expression for temp (computed from variables).
    if (!cir.elaborate({}, "__topdef__", "__topinst__", { PTParameters(PVv( PV("reltol", 1e-4))), p.parseParameters("temp=var1") }, nullptr, s)) {
        Simulator::err() << "Elaboration failed.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    // Print circuit option temp
    std::cout << "temp=" << cir.simulatorOptions().core().temp << "\n\n";
    
    cir.dumpHierarchy(0, Simulator::out());

    // Analysis description
    auto dcDesc = PTAnalysis("dc1", "op");
    dcDesc
        // Outermost sweep first
        .add(PTSweep("mode")
            .add(PV("instance", "x1"))
            .add(PV("parameter", "mode"))
            .add(PV("values", IntVector{0, 1}))
        );

    // Analysis object, no options map (options that are given as parameterized expressions)
    auto dc = Analysis::create(dcDesc, cir, s);
    if (!dc) {
        Simulator::err() << "Failed to create analysis.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    // Run analysis
    auto [ok, canResume] = dc->run(s);
    if (!ok) {
        Simulator::err() << "Analysis failed.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }
    Simulator::out() << "Analysis OK. Can resume: " << (canResume ? "true" : "false") << "\n";

    // Cleanup
    delete dc;

    // Run postprocessing
    runProcess(pythonBinary, {"runme.py"}, &pythonLibraryPath, false, false);

    return 0;
}
