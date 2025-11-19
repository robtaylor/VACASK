#include "libplatform.h"
#include "simulator.h"
#include "parser.h"
#include "openvafcomp.h"
#include "circuit.h"
#include "processutils.h"

using namespace sim;

// Demonstrate elaboration that changes the topology
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
    ParserTables tab("Variables and parameterized options sweep");

    // Parser, needs tables to store stats when parsing expressions and parameters
    Parser p(tab);

    // Title, loads, default ground (0)
    tab
        .add(PTLoad("resistor.osdi"))
        .add(PTLoad("diode.osdi"))
        .defaultGround()
        .setDefaultSubDef(
            // Toplevel subcircuit definition, no name, no terminals
            PTSubcircuitDefinition()
            // Resistor model res, capacitor model cap, voltage source model vsrc
            .add(PTModel("res", "resistor"))
            .add(PTModel("dio", "diode")
                .add(p.parseParameters("is=1e-12 n=2 rs=1 eg=1.2 xti=2"))
            )
            .add(PTModel("vsrc", "vsource"))
            // First topology: R-D, no terminals (toplevel circuit)
            .add(PTSubcircuitDefinition("sub1")
                .add(PTInstance("r1", "res", {"1", "2"}).add(PV("r", 10)))
                .add(PTInstance("d1", "dio", {"2", "0"}))
            )
            // Second topology: D-R, no terminals (toplevel circuit)
            .add(PTSubcircuitDefinition("sub2")
                .add(PTInstance("d1", "dio", {"1", "2"}))
                .add(PTInstance("r1", "res", {"2", "0"}).add(PV("r", 10)))
            )
            // Power supply in default toplevel circuit
            .add(PTInstance("v1", "vsrc", {"1", "0"}).add(PV("dc", 10)))
        )
        // Embedded Python postprocessing script
        .add(PTEmbed("runme.py", R"script(
from rawfile import rawread
import numpy as np
import matplotlib.pyplot as plt 
import sys

dc1 = rawread('dc1.raw').get()
dc2 = rawread('dc2.raw').get()
print("dc1 v(2)=", dc1["2"])
print("dc2 v(2)=", dc2["2"])
print("Sum (should be 10):", dc1["2"]+dc2["2"])
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

    // Elaborate default toplevel circuit + sub1
    // If not specified otherwise before first elaboration, 
    // default simulator options are used. 
    if (!cir.elaborate({"sub1"}, "__topdef__", "__topinst__", nullptr, s)) {
        Simulator::err() << "Elaboration failed.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    // Analysis description (sweep temp via variable)
    auto dc1Desc = PTAnalysis("dc1", "op");
    
    // Analysis object, no saves, with an options map
    auto dc1 = Analysis::create(dc1Desc, cir, s);
    if (!dc1) {
        Simulator::err() << "Failed to create analysis.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    // Run analysis
    if (auto [ok, canResume] = dc1->run(s); !ok) {
        Simulator::err() << "DC1 analysis failed.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    } else {
        Simulator::out() << "DC1 analysis OK. Can resume: " << (canResume ? "true" : "false") << "\n";
    }

    // Elaboration does not reset variables and options. 
    // Need to reset them manually by calling clearVariables() and clearOptions(). 
    // Elaborate default toplevel circuit + sub2
    if (!cir.elaborate({"sub2"}, "__topdef__", "__topinst__", nullptr, s)) {
        Simulator::err() << "Elaboration failed.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    // Analysis description (sweep temp via variable)
    auto dc2Desc = PTAnalysis("dc2", "op");
    
    // Analysis object, no saves, with an options map
    auto dc2 = Analysis::create(dc2Desc, cir, s);
    if (!dc2) {
        Simulator::err() << "Failed to create analysis.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    // Run analysis
    if (auto [ok, canResume] = dc2->run(s); !ok) {
        Simulator::err() << "DC2 analysis failed.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    } else {
        Simulator::out() << "DC2 analysis OK. Can resume: " << (canResume ? "true" : "false") << "\n";
    }

    // Cleanup
    delete dc1, dc2;

    // Run postprocessing
    runProcess(pythonBinary, {"runme.py"}, &pythonLibraryPath, false, false);

    return 0;
}

// TODO:
//   expression for analysis/sweep parameter, sweep it
//   change parameter, elaborate changes
//   change option, elaborate changes
//   change variable that affects an instance, elaborate changes
