#include "libplatform.h"
#include "simulator.h"
#include "parser.h"
#include "openvafcomp.h"
#include "circuit.h"
#include "processutils.h"

using namespace sim;

// Demonstrate partial elaboration
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
            .add(PTModel("res", "resistor"))
            .add(PTModel("dio", "diode")
                .add(p.parseParameters("is=1e-12 n=2 rs=1 eg=1.2 xti=2"))
            )
            .add(PTModel("vsrc", "vsource"))

            .add(PTBlockSequence()
                // Fixed resistor
                .add(p.parseExpression("fixed!=0"), PTBlock()
                    .add(PTInstance("r1", "res", {"1", "2"})
                        .add(PV("r", 100))
                    )
                )
                // Else block has an empty Rpn object as condition
                .add(Rpn(), PTBlock()
                    .add(PTInstance("r1", "res", {"1", "2"})
                        // A formula is used for computing resistance that depends on circuit temperature, TC=1%/degC
                        // Nominal temperature of the resistor is set by the tnominal circuit variable. 
                        .add(p.parseParameters("r=10*(1+($temp-tnominal)/100)"))
                    )
                )
            )
            .add(PTInstance("d1", "dio", {"2", "0"}))
            .add(PTInstance("v1", "vsrc", {"1", "0"}).add(PV("dc", 10)))
        )
        // Embedded Python postprocessing script
        .add(PTEmbed("runme.py", R"script(
from rawfile import rawread
import numpy as np
import matplotlib.pyplot as plt 
import sys

dc1 = rawread('dc1.raw').get()
fig1, ax1 = plt.subplots(1, 1, figsize=(6,4), dpi=100, constrained_layout=True)
fig1.axes[0].plot(dc1["tnom"], dc1["2"], "r", marker=".")

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

    // Lambda fopr printing the state
    auto state = [&cir]() {
        Simulator::out() << "reltol=" << cir.simulatorOptions().core().reltol << "\n";
        auto var1 = cir.getVariable("tnominal");
        Simulator::out() << "tnominal=" << (var1?var1->str():"not found") << "\n";
        auto var2 = cir.getVariable("fixed");
        Simulator::out() << "fixed=" << (var2?var2->str():"not found") << "\n";
        auto [ok1, val1] = cir.instanceParameter("v1", "dc");
        Simulator::out() << "v1 dc=" << (ok1?val1.str():"not found") << "\n";
        auto [ok2, val2] = cir.instanceParameter("r1", "r");
        Simulator::out() << "r1 r=" << (ok2?val2.str():"not found") << "\n";
    };

    // Lambda for partial elaboration
    auto partialElaborate = [&cir]() {
        Simulator::out() << "\nPartial elaboration\n";

        Status s;
        if (auto [ok, topologyChange, bindingNeeded] = cir.elaborateChanges(nullptr, s); !ok) {
            Simulator::err() << "Partial elaboration failed.\n";
            Simulator::err() << s.message() << "\n";
            exit(1);
        } else {
            // Did the topology change? 
            if (topologyChange) {
                Simulator::out() << "Topology changed.\n";
            }
            // Binding is performed during a sweep inside the analysis. 
            // No need to do anything after partial elaboration at API level. 
            if (bindingNeeded) {
                Simulator::out() << "Analysis rebinding needed.\n";
            }
        }
    };

    // Set circuit variables (needed at elaboration)
    cir.setVariable("tnominal", 27);
    cir.setVariable("fixed", 1);

    Simulator::out() << "\nBefore elaboration\n";
    state();

    // Elaborate default toplevel circuit
    if (!cir.elaborate({}, "__topdef__", "__topinst__", nullptr, s)) {
        Simulator::err() << "Elaboration failed.\n";
        Simulator::err() << s.message() << "\n";
        exit(1);
    }
    Simulator::out() << "\nAfter elaboration\n";
    state();

    // Change some circuit parameters, variables, and hierarchical parameters
    cir.setOption("reltol", 1e-6);
    cir.setVariable("tnominal", 40);
    cir.setInstanceParameter("v1", "dc", 20);

    // Elaborate changes
    partialElaborate();
    Simulator::out() << "After first partial elaboration\n";
    state();
    
    // Switch to temperature-dependent resistor
    cir.setVariable("fixed", 0);
    
    // Elaborate changes
    // After elaboration the whole affected subcircuit is rebuilt. 
    // Therefore v1 dc is reset to its netlist value. 
    partialElaborate();
    Simulator::out() << "After second partial elaboration\n";
    state();

    // Set nominal temperature to 27
    cir.setVariable("tnominal", 27);
    // Resistor should now be r=10
    partialElaborate();
    Simulator::out() << "After third partial elaboration\n";
    state();

    // Analysis description (sweep tnominal variable)
    auto dc1Desc = PTAnalysis("dc1", "op");
    dc1Desc
        .add(PTSweep("tnom")
            .add(PV("variable", "tnominal"))
            .add(PV("from", 0))
            .add(PV("to", 100))
            .add(PV("step", 2))
        );
    
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

    // Because R is above D and increasing tnominal decreases the resistance
    // the output voltage should increase when tnominal increases. 

    // Cleanup
    delete dc1;

    // Run postprocessing
    runProcess(pythonBinary, {"runme.py"}, &pythonLibraryPath, false, false);

    return 0;
}
