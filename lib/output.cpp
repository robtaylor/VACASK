#include "output.h"
#include "common.h"
#include <limits>
#include <complex>


namespace NAMESPACE {

using namespace std::complex_literals;

Complex OutputSource::defaultBad = 0.0 + 0.0i; 

OutputDescriptor::OutputDescriptor(Type t, Id outputName) : type(t), name(outputName) {
}

OutputDescriptor::OutputDescriptor(Type t, Id outputName, size_t n) : type(t), name(outputName)  {
    if (n>std::numeric_limits<NdxType>::max()) {
        throw std::length_error("Output descriptor, index too large.");
    }
    ndx = n;
}

OutputDescriptor::OutputDescriptor(Type t, Id outputName, Id id) : type(t), name(outputName), id(id)  {
}

OutputDescriptor::OutputDescriptor(Type t, Id outputName, Id id1, Id id2) : type(t), name(outputName)  {
    idId.id1 = id1;
    idId.id2 = id2;
}

OutputDescriptor::OutputDescriptor(Type t, Id outputName, Id id, size_t n) : type(t), name(outputName)  {
    if (n>std::numeric_limits<NdxType>::max()) {
        throw std::length_error("Output descriptor, index too large.");
    }
    idNdx.id = id;
    idNdx.ndx = n;
}


/*
Id OutputTdRhs::name() const {
    return node_;
}

void OutputTdRhs::dump(std::ostream& os) {
    os << "RHStd("+std::string(node_)+")";
}


Id OutputOutvar::name() const {
    return Id(std::string(instance_)+"."+std::string(outvar_));
}

void OutputOutvar::dump(std::ostream& os) {
    os << "Outvar("+std::string(instance_)+","+std::string(outvar_)+")";
}


Id OutputSweepvar::name() const {
    return name_;
}

Sweeper::SweepIndexType OutputSweepvar::index() const {
    return ndx;
}

void OutputSweepvar::dump(std::ostream& os) {
    os << "Sweepvar(" << name_ << ")";
}
*/

Output::Output(const std::string& baseName, DescriptorList& descriptors, SourcesList& sources)
    : baseName_(baseName), descrs(descriptors), srcs(sources), 
    title_("Output file"), plotname_("analysis") {
}

}
