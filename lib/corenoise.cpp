#include <iomanip>
#include <cmath>
#include <filesystem>
#include "corenoise.h"
#include "simulator.h"
#include "answeep.h"
#include "context.h"
#include "common.h"
#include <numbers>

namespace NAMESPACE {

// Default parameters
NoiseParameters::NoiseParameters() {
    opParams.writeOutput = 0;
}

template<> int Introspection<NoiseParameters>::setup() {
    registerMember(out);
    registerMember(in);
    registerMember(from);
    registerMember(to);
    registerMember(step);
    registerMember(mode);
    registerMember(points);
    registerMember(values);
    registerMember(dumpop);
    registerNamedMember(opParams.nodeset, "nodeset");
    registerNamedMember(opParams.store, "store");

    return 0;
}
instantiateIntrospection(NoiseParameters);


NoiseCore::NoiseCore(
    OutputDescriptorResolver& parentResolver, NoiseParameters& params, OperatingPointCore& opCore, 
    std::unordered_map<std::pair<Id, Id>, size_t>& contributionOffset, 
    Circuit& circuit, 
    KluRealMatrix& dcJacobian, VectorRepository<double>& dcSolution, VectorRepository<double>& dcStates, 
    KluComplexMatrix& acMatrix, Vector<Complex>& acSolution, 
    Vector<double>& results, double& powerGain, double& outputNoise
) : AnalysisCore(parentResolver, circuit), params(params), outfile(nullptr), opCore_(opCore), 
    dcSolution(dcSolution), dcStates(dcStates), dcJacobian(dcJacobian), 
    acMatrix(acMatrix), acSolution(acSolution), 
    contributionOffset(contributionOffset), 
    results(results), powerGain(powerGain), outputNoise(outputNoise) {
    
    // Set analysis type for the initial operating point analysis
    auto& elsSystem = opCore_.solver().evalSetup();
    elsSystem.staticAnalysis = true;
    elsSystem.dcAnalysis = false;
    elsSystem.noiseAnalysis = true;
}

NoiseCore::~NoiseCore() {
    delete outfile;
}

bool NoiseCore::resolveOutputDescriptors(bool strict) {
    clearError();
    // Clear contribution offsets
    contributionOffset.clear();
    size_t atOffset = 2;
    // Clear results
    results.clear();
    // Resolve output descriptors
    bool ok = true; 
    for (auto it = outputDescriptors.cbegin(); it != outputDescriptors.cend(); ++it) {
        Id name;
        Id contrib;
        ParameterIndex contribIndex;
        bool found;
        size_t ndx;
        Instance *inst;
        switch (it->type) {
            case OutdNoiseContribInst:
                name = it->id;
                // Find instance
                inst = circuit.findInstance(name);
                if (strict) {
                    if (!inst) {
                        setError(NoiseError::NotFound);
                        errorInstance = name;
                        ok = false;
                        break;
                    }
                }
                if (inst) {
                    // Add to contribution offsets
                    auto [insIt, inserted] = contributionOffset.insert({{name, Id()}, contributionOffset.size()});
                    outputSources.emplace_back(&results, insIt->second);
                } else {
                    // Instance not found, constant source
                    outputSources.emplace_back();
                }
                break;
            case OutdNoiseContribInstPartial:
                name = it->idId.id1;
                contrib = it->idId.id2;
                // Find instance
                inst = circuit.findInstance(name);
                if (strict && !inst) {
                    setError(NoiseError::NotFound);
                    errorInstance = name;
                    ok = false;
                    break;
                }
                // Find contrib
                if (inst) {
                    std::tie(contribIndex, found) = inst->noiseSourceIndex(contrib);
                    if (strict && !found) {
                        setError(NoiseError::ContribNotFound);
                        errorInstance = name;
                        errorContrib = contrib;
                        ok = false;
                        break;
                    }
                    if (found) {
                        // Add to contribution offsets
                        auto [insIt, inserted] = contributionOffset.insert({{name, contrib}, contributionOffset.size()});
                        outputSources.emplace_back(&results, insIt->second);
                    } else {
                        // Contribution not found, constant source
                        outputSources.emplace_back();
                    }
                } else {
                    // Instance not found, constant source
                    outputSources.emplace_back();
                }
                break;
            case OutdFrequency:
                outputSources.emplace_back(&frequency);
                break;
            case OutdOutputNoise:
                outputSources.emplace_back(&outputNoise);
                break;
            case OutdPowerGain:
                outputSources.emplace_back(&powerGain);
                break;
            default:
                // Delegate to parent
                ok = parentResolver.resolveOutputDescriptor(*it, outputSources, strict);
        }
        if (!ok) {
            break;
        }
    }
    return ok;
}

bool NoiseCore::addCoreOutputDescriptors() {
    clearError();
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    if (!addOutputDescriptor(OutputDescriptor(OutdFrequency, "frequency"))) {
        lastError = Error::Descriptor;
        errorId = "frequency";
        return false;
    }
    if (!addOutputDescriptor(OutputDescriptor(OutdOutputNoise, "onoise"))) {
        lastError = Error::Descriptor;
        errorId = "onoise";
        return false;
    }
    if (!addOutputDescriptor(OutputDescriptor(OutdPowerGain, "gain"))) {
        lastError = Error::Descriptor;
        errorId = "gain";
        return false;
    }
    return true;
}

bool NoiseCore::addDefaultOutputDescriptors() {
    // If output is suppressed, skip all this work
    if (!params.writeOutput) {
        return true;
    }
    if (savesCount==0) {
        // Add total noise contributions of all instances (details=false)
        return addAllNoiseContribInst(PTSave(Loc::bad, "default", Id(), Id()), false);
    }
    return true;
}

bool NoiseCore::initializeOutputs(Id name, Status& s) {
    // Create output file if not created yet
    if (!outfile) {
        outfile = new OutputRawfile(
            name, outputDescriptors, outputSources,
            (circuit.simulatorOptions().core().rawfile==SimulatorOptions::rawfileBinary ? OutputRawfile::Flags::Binary : OutputRawfile::Flags::None) |
                OutputRawfile::Flags::Padded);
        outfile->setTitle(circuit.title());
        outfile->setPlotname("Small-Signal Noise Analysis");
    }
    outfile->prologue();

    return true;
}

bool NoiseCore::finalizeOutputs(Status& s) {
    outfile->epilogue();
    delete outfile;
    outfile = nullptr;
    return true;
}

bool NoiseCore::deleteOutputs(Id name, Status& s) {
    if (!params.writeOutput) {
        return true;
    }

    // Cannot assume outfile is available
    auto fname = std::string(name)+".raw";
    if (std::filesystem::exists(fname)) {
        std::filesystem::remove(fname);
    }
    return true;
}
    
bool NoiseCore::rebuild(Status& s) {
    // AC analysis matrix
    if (!acMatrix.rebuild(circuit.sparsityMap(), circuit.unknownCount())) {
        acMatrix.formatError(s);
        return false;
    }
    
    // Resistive Jacobian entries remain bound to OP Jacobian, 
    // reactive parts will be bound to imaginary entries of acMatrix
    if (!circuit.bind(nullptr, Component::Real, std::nullopt, &acMatrix, Component::Imaginary, std::nullopt, s)) {
        return false;
    }
    
    return true;
}

// System of equations is 
//   (G(x) + i C(x)) dx = dJ
CoreCoroutine NoiseCore::coroutine(bool continuePrevious) {
    acMatrix.setAccounting(circuit.tables().accounting());
    
    clearError();
    auto n = circuit.unknownCount(); 
    // Make sure structures are large enough
    acSolution.resize(n+1);
    results.resize(contributionOffset.size());
    zero(results);
    
    // Get output unknowns
    auto [outOk, up, un] = getOutput(params.out);
    if (!outOk) {
        co_yield CoreState::Aborted;
    }

    // Get input source
    auto [inOk, inputSource] = getInput(params.in);
    if (!inOk) {
        co_yield CoreState::Aborted;
    }
    
    // Compute operating point
    auto opOk = opCore_.run(continuePrevious);
    if (!opOk) {
        setError(NoiseError::OperatingPointError);
        co_yield CoreState::Aborted;
    }

    auto& options = circuit.simulatorOptions().core();
    Int debug = options.smsig_debug;

    if (debug>0) {
        Simulator::dbg() << "Starting small-signal noise analysis.\n";
    }
    
    // Evaluate resistive and reactive Jacobian, evaluate noise, bypass is not allowed
    EvalSetup esReactNoise { 
        // Inputs, can be set here (we do not rotate)
        .solution = &dcSolution, 
        .states = &dcStates, 
        
        // Evaluation type reported to the model
        .noiseAnalysis = true, 

        // Evaluation 
        .enableLimiting = false, 
        .evaluateResistiveJacobian = true, 
        .evaluateReactiveJacobian = true, 
        .evaluateNoise = true, 
    };

    LoadSetup lsReacNoise { 
        // Outputs
        .loadResistiveJacobian = false, 
        .loadReactiveJacobian = true, 
    };

    // Copy OP Jacobian to real part of acMatrix, zero out imaginary part. 
    // Because the real part is taken from OP Jacobian it includes
    // the shunt resistors. 
    auto nnz = dcJacobian.nnz();
    auto Jr = dcJacobian.data();
    auto M = acMatrix.data();
    for(decltype(nnz) i=0; i<nnz; i++) {
        M[i] = Jr[i];
    }
    
    // Evaluate Jacobians 
    // Actually we only need to evaluate the reactive Jacobian 
    // because the resistive part was evaluated by OP analysis
    // We do both here in case OpenVAF has bugs with this corner case :)
    if (!circuit.evalAndLoad(&esReactNoise, nullptr, nullptr)) {
        // Load error
        setError(NoiseError::EvalAndLoad);
        if (debug>0) {
            Simulator::dbg() << "Error in AC Jacobian / noise evaluation.\n";
        }
        co_yield CoreState::Aborted;
    }

    // Handle Abort, Finish, Stop
    if (esReactNoise.requests.abort) {
        if (debug>0) {
            Simulator::dbg() << "Abort requested during AC Jacobian / noise  evaluation. Exiting.\n";
        }
        co_yield CoreState::Aborted;
    }
    if (esReactNoise.requests.finish) {
        if (debug>0) {
            Simulator::dbg() << "Finish requested during AC Jacobian / noise evaluation. Exiting.\n";
        }
        co_yield CoreState::Finished;
    }
    if (esReactNoise.requests.stop) {
        if (debug>0) {
            Simulator::dbg() << "Stop requested during AC Jacobian / noise evaluation. Exiting.\n";
        }
        co_yield CoreState::Stopped;
    }

    // Create sweeper, put it in unique ptr to free it when method returns
    ScalarSweep sweeper;
    if (!sweeper.setup(params, errorStatus)) {
        setError(NoiseError::Sweeper);
        co_yield CoreState::Aborted;
    }
    if (progressReporter) {
        progressReporter->setValueFormat(ProgressReporter::ValueFormat::Scientific, 6);
        progressReporter->setValueDecoration("", "Hz");    
    }
    initProgress(sweeper.count(), 0);
    
    // Frequency sweep
    sweeper.reset();
    bool finished = false;
    frequency = -1.0;
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4);
    bool error = false;
    do {
        // Compute should always succeed
        Value v;
        if (!sweeper.compute(v, errorStatus)) {
            setError(NoiseError::SweepCompute);
            error = true;
            break;
        }

        // The value, however, must be convertible to real
        if (!v.convertInPlace(Value::Type::Real, errorStatus)) {
            setError(NoiseError::SweepCompute);
            error = true;
            break;
        }
        frequency = v.val<Real>();
        double omega = 2*std::numbers::pi*frequency;
        
        if (debug>0) {
            ss.str(""); ss << frequency;
            Simulator::dbg() << "frequency=" << ss.str() << "\n";
        }

        // Load AC matrix, we must update the imaginary part only
        acMatrix.zero(Component::Imaginary);
        lsReacNoise.reactiveJacobianFactor = omega;
        if (!circuit.evalAndLoad(nullptr, &lsReacNoise, nullptr)) {
            // Load error
            setError(NoiseError::EvalAndLoad);
            if (debug>0) {
                Simulator::dbg() << "Error in AC Jacobian load.\n";
            }
            error = true;
            break;
        }

        if (debug>=101) {
            Simulator::dbg() << "Linear system matrix\n";
            acMatrix.dump(Simulator::dbg()); 
            Simulator::dbg() << "\n";
        }

        // Check if matrix entries are finite, no need to check RHS 
        // since we loaded it without any computation (i.e. we only used mag and phase)
        if (options.matrixcheck && !acMatrix.isFinite(true, true)) {
            setError(NoiseError::MatrixError);
            if (debug>2) {
                Simulator::dbg() << "A matrix entry is not finite.\n";
            }
            error = true;
            break;
        }

        // Factor
        bool forceFullFactorization = false;        
        if (acMatrix.isFactored()) {
            // Refactor (if possible)
            if (!acMatrix.refactor()) {
                // Failed, try again by fully factoring
                forceFullFactorization = true;
            } 
        }
        if (forceFullFactorization || !acMatrix.isFactored()) {
            // Full factorization
            if (!acMatrix.factor()) {
                // Failed, give up
                setError(NoiseError::MatrixError);
                if (debug>0) {
                    Simulator::dbg() << "LU factorization failed.\n";
                }
                error = true;
                break;
            }
        }
        // Check if matrix is singular
        double rcond;
        if (!acMatrix.rcond(rcond)) {
            setError(NoiseError::MatrixError);
            if (debug>0) {
                Simulator::dbg() << "Condition number estimation failed.\n";
            }
            error = true;
            break;
        }
        if (rcond<1e-15) {
            if (debug>0) {
                Simulator::dbg() << "Matrix is close to singular.\n";
            }
            setError(NoiseError::SingularMatrix);
            error = true;
            break;
        }

        
        // Compute power gain
        zero(acSolution); 
        auto [e1, e2] = inputSource->sourceExcitation(circuit);
        acSolution[e1] += -inputSource->scaledUnityExcitation();
        acSolution[e2] -= -inputSource->scaledUnityExcitation();

        if (debug>=100) {
            Simulator::dbg() << "Linear system for power gain\n";
            acMatrix.dump(Simulator::dbg(), dataWithoutBucket(acSolution)); 
            Simulator::dbg() << "\n";
        }

        // Solve, set bucket to 0.0
        if (!acMatrix.solve(dataWithoutBucket(acSolution))) {
            setError(NoiseError::MatrixError);
            if (debug>2) {
                Simulator::dbg() << "Failed to solve factored system.\n";
            }
            error = true;
            break;
        }
        acSolution[0] = 0.0;

        if (options.solutioncheck && !acMatrix.isFinite(dataWithoutBucket(acSolution), true, true)) {
            setError(NoiseError::SolutionError);
            if (options.smsig_debug) {
                Simulator::dbg() << "A solution entry is not finite. Solver failed.\n";
            }
            error = true;
            break;
        }
        
        // Power gain
        // Note that for $mfactor!=1 this gain is from the magnitude of the source 
        // to the given output, not from total source value to the output. 
        // $mfactor does not change anything for a gain computed from a voltage source, 
        // but it affects the gain computed from a current source. 
        auto tf = (acSolution[up] - acSolution[un]);
        powerGain = std::abs(tf);
        powerGain *= powerGain;

        Vector<double> noiseDensity;
        
        // Set total output noise to 0
        outputNoise = 0.0;

        // Go through all instances
        auto ndev = circuit.deviceCount();
        for(decltype(ndev) idev=0; idev<ndev; idev++) {
            auto dev = circuit.device(idev);
            auto nmod = dev->modelCount();
            for(decltype(nmod) imod=0; imod<nmod; imod++) {
                auto mod = dev->model(imod);
                auto ninst = mod->instanceCount();
                for(decltype(ninst) iinst=0; iinst<ninst; iinst++) {
                    auto inst = mod->instance(iinst);

                    // Skip instances without noise sources
                    auto nSources = inst->noiseSourceCount(); 
                    if (nSources<=0) {
                        continue;
                    }

                    // Instance name
                    auto name = inst->name();

                    if (debug>1) {
                        Simulator::dbg() << "  instance '" << std::string(name) << "'\n";
                    }

                    // Collect noise excitations
                    noiseDensity.resize(nSources);
                    if (!inst->loadNoise(circuit, frequency, noiseDensity.data())) {
                        setError(NoiseError::PsdError);
                        if (debug>0) {
                            Simulator::dbg() << "Failed to compute noise.\n";
                        }
                        error = true;
                        break;
                    }

                    // Go through all noise sources
                    double sourceContribution = 0.0;
                    double totalInstanceContribution = 0.0;
                    for(decltype(nSources) ndx=0; ndx<nSources; ndx++) {
                        auto contrib = inst->noiseSourceName(ndx);

                        if (debug>1) {
                            Simulator::dbg() << "    contribution '" << std::string(contrib) << "'\n";
                        }

                        // Compute gain from noise source to output
                        zero(acSolution); 
                        auto [e1, e2] = inst->noiseExcitation(circuit, ndx);

                        // Set RHS, load negated unity excitation to get the true value of response after solve()
                        // Here this is not neccessary because we are working with the absolute value of the response. 
                        acSolution[e1] += -1.0;
                        acSolution[e2] -= -1.0;

                        if (debug>=100) {
                            Simulator::dbg() << "Linear system for contribution '"+std::string(contrib)+"' of '"+std::string(name)+"'\n";
                            acMatrix.dump(Simulator::dbg(), dataWithoutBucket(acSolution)); 
                            Simulator::dbg() << "\n";
                        }

                        // Solve, set bucket to 0.0
                        if (!acMatrix.solve(dataWithoutBucket(acSolution))) {
                            setError(NoiseError::MatrixError);
                            if (debug>2) {
                                Simulator::dbg() << "Failed to solve factored system.\n";
                            }
                            error = true;
                            break;
                        }
                        acSolution[0] = 0.0;
                        
                        // Power gain from noise source to output
                        auto tf = acSolution[up] - acSolution[un];
                        auto gain = std::abs(tf);
                        gain *= gain;
                        
                        // Contribution - take absolute value of noise power spectral density (it may be negative). 2
                        // See Coram et. al., Flicker Noise Formulations in Compact Models, IEEE TCAD, vol 39, 2020. 
                        sourceContribution = gain * std::abs(noiseDensity[ndx]);
                        totalInstanceContribution += sourceContribution;
                        // Simulator::dbg() << "       src=" << sourceContribution << "  inst=" << totalInstanceContribution << "\n";

                        // std::cout << "gain = " << gain << "\n";
                        // std::cout << "src  = " << noiseDensity[ndx] << "\n";
                        // std::cout << "out  = " << instanceContribution << "\n";

                        // Store instance contribution
                        auto it = contributionOffset.find({name, contrib});
                        if (it!=contributionOffset.end()) {
                            // std::cout << "store: " << it->second << "\n";
                            results[it->second] = sourceContribution;
                        }
                    }
                    // End of noise sources loop

                    // Store total instance contribution
                    auto it = contributionOffset.find({name, Id()});
                    if (it!=contributionOffset.end()) {
                        // std::cout << "store: " << it->second << "\n";
                        results[it->second] = totalInstanceContribution;
                    }

                    // Add to output noise
                    outputNoise += totalInstanceContribution;
                    // std::cout << "total = " << outputNoise << "\n";

                    if (error) {
                        break;
                    }
                }
                // End of instances loop

                if (error) {
                    break;
                }
            }
            // End of models loop

            if (error) {
                break;
            }
        }
        // End of devices loop
        
        if (error) {
            break;
        }

        // std::cout << "total = " << outputNoise << "\n";

        // Dump solution
        if (params.writeOutput && outfile) {
            outfile->addPoint();
        }

        finished = sweeper.advance();

        setProgress(sweeper.at(), frequency);
    } while (!finished && !error);

    if (debug>0) {
        Simulator::dbg() << "Noise frequency sweep " << (finished ? "completed" : "exited prematurely") << ".\n";
    }

    if (!finished) {
        errorFreq = frequency;
    }

    // No need to bind resistive Jacobian entries. 
    // OP analysis will still work fine, even in sweep. 
    // We only changed the bindings of the reactive Jacobian entries. 
    
    if (finished) {
        co_yield CoreState::Finished;
    } else {
        co_yield CoreState::Aborted;
    }
}

bool NoiseCore::run(bool continuePrevious) {
    auto c = coroutine(continuePrevious);
    bool ok = true;
    while (!c.done()) {
        if (c.resume()==CoreState::Aborted) {
            ok = false;
            break;
        };
    }
    return ok;
}

bool NoiseCore::formatError(Status& s) const {
    auto nr = UnknownNameResolver(circuit);
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4);
    
    // First, handle AnalysisCore errors
    if (lastError!=Error::OK) {
        AnalysisCore::formatError(s);
        return false;
    }
    
    // Then handle NoiseCore errors
    switch (lastNoiseError) {
        case NoiseError::NotFound:
            s.set(Status::Analysis, std::string("Instance '")+std::string(errorInstance)+"' not found.");
            break;
        case NoiseError::ContribNotFound:
            s.set(Status::Analysis, std::string("Noise contribution '")+std::string(errorContrib)+"' of instance '"+std::string(errorInstance)+"' not found.");
            break;
        case NoiseError::Sweeper:
        case NoiseError::SweepCompute:
            s.set(errorStatus);
            break;
        case NoiseError::EvalAndLoad:
            s.set(Status::Analysis, "Jacobian evaluation failed.");
            break;
        case NoiseError::PsdError:
            s.set(Status::Analysis, "Power spectral density evaluation failed.");
            break;
        case NoiseError::MatrixError:
            acMatrix.formatError(s, &nr);
            break;
        case NoiseError::SolutionError:
            acMatrix.formatError(s, &nr);
            s.extend("Solution component is not finite.");
            break;
        case NoiseError::OperatingPointError:
            opCore_.formatError(s);
            break;
        case NoiseError::SingularMatrix:
            s.set(Status::Analysis, "Matrix is close to singular.");
            break;
        case NoiseError::BadFrequency:
            s.set(Status::Analysis, "Frequency value cannot be converted to real.");
            break;
        default:
            return true;
    }
    if (errorFreq>=0) {
        ss.str(""); ss << errorFreq;
        s.extend(std::string("Leaving frequency sweep at frequency=")+ss.str()+".");
    } else {
        s.extend("Leaving frequency sweep.");
    }
    return false;
}


void NoiseCore::dump(std::ostream& os) const {
    AnalysisCore::dump(os);
    os << "  Results\n";
    os << "    Output noise: " << outputNoise << "\n";
    os << "    Power gain: " << powerGain << "\n";
    for(auto& it : contributionOffset) {
        auto [inst, contrib] = it.first;
        auto ndx = it.second;
        if (contrib) {
            os << "    n("+std::string(inst)+","+std::string(contrib)+") " << results[ndx] << "\n";
        } else {
            os << "    n("+std::string(inst)+") " << results[ndx] << "\n";
        }
    }
}

}
