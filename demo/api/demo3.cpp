#include "libplatform.h"
#include "simulator.h"
#include "parser.h"
#include "openvafcomp.h"
#include "circuit.h"
#include "processutils.h"

using namespace sim;

// Demonstrate circuit variables and parameterized options
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
            // Element r1, master is res, terminals 1 and 2, constant parameter r=1000, no parameter expressions
            .add(PTInstance("r1", "res", {"1", "2"})
                .add(PV{"r", 100})
            )
            .add(PTInstance("d1", "dio", {"2", "0"}))
            // Pulse source
            .add(PTInstance("v1", "vsrc", {"1", "0"})
                // Add parsed parameters
                .add(PV("dc", 2))
            )
        )
        // Embedded Python postprocessing script
        .add(PTEmbed("runme.py", R"script(
from rawfile import rawread
import numpy as np
import matplotlib.pyplot as plt 
import sys

dc1 = rawread('dc1.raw').get()
dc2 = rawread('dc2.raw').get()
print("Vectors:", dc1.names)
fig1, ax1 = plt.subplots(1, 1, figsize=(6,4), dpi=100, constrained_layout=True)
fig1.axes[0].plot(dc2["temp"], dc2["2"], "g")
fig1.axes[0].plot(dc1["temp"], dc1["2"], "r", marker=".", linestyle="none")

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

    // Define myvar circuit variable
    // Note that elaboration does not clear circuit variables 
    // because variables are used for evaluating options expressions 
    // at elaboration. 
    // If you want to clear the variables, call clearVariables() manually. 
    cir.setVariable("myvar", 0);

    // Elaborate it, just the default toplevel subcircuit 
    // (empty list of subcircuit definition names). 
    // Use __topdef__ and __topinst__ as prefixes for toplevel 
    // subcircuit model and toplevel subcircuit instance names. 
    // Do not collect device requests (i.e. abort/finish/stop). 
    cir.setOption("reltol", 1e-4);
    if (!cir.elaborate({}, "__topdef__", "__topinst__", nullptr, s)) {
        Simulator::err() << "Elaboration failed.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }

    cir.dumpHierarchy(0, Simulator::out());


    // Analysis description (sweep temp via variable)
    auto dc1Desc = PTAnalysis("dc1", "op");
    dc1Desc
        .add(PTSweep("temp")
            .add(PV("variable", "myvar"))
            .add(PV("from", -50))
            .add(PV("to", 100))
            .add(PV("step", 2))
        );
    
    // Options map that is applied in analysis
    // Holds only options that are defined as expressions. 
    
    // parseParameters() creates a PTParameters which holds actual vaues and expressions. 
    auto anPar = p.parseParameters(
        // temp will be computed from an expression that depends on the myvar circuit variable
        // Any option specified as a constant will be ignored. 
        "temp=myvar"
    );
    
    // Analysis object, no saves, with parameterized options
    auto dc1 = Analysis::create(dc1Desc, cir, s);
    // Add an analysis parameter expression. 
    // This one will be affected because we sweep myvar.  
    dc1->add(anPar);
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

    // After the analysis is finished the circuit state 
    // (instance and model parameters, options, variables, topology) 
    // should be the same as before analysis. 
    
    // Analysis description (sweep temp directly)
    auto dc2Desc = PTAnalysis("dc2", "op");
    dc2Desc
        .add(PTSweep("temp")
            .add(PV("option", "temp"))
            .add(PV("from", -50))
            .add(PV("to", 100))
            .add(PV("step", 2))
        );
    
    // Analysis object, no saves, no options map
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

    // After the analysis is finished the circuit state 
    // (instance and model parameters, options, variables, topology) 
    // should be the same as before analysis. 
    
    // Cleanup
    delete dc1, dc2;

    // Run postprocessing
    runProcess(pythonBinary, {"runme.py"}, &pythonLibraryPath, false, false);

    return 0;
}
