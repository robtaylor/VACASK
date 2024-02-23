#include "node.h"

namespace NAMESPACE {

Node::Node(Id name, Flags f) 
    : name_(name), FlagBase<NodeFlags>(f), refCnt_(0) {
}

Node::~Node() {
}

RefCountIndex Node::incRef() {
    if (refCnt_==std::numeric_limits<decltype(refCnt_)>::max()) {
        throw std::length_error("Too many references to a node.");
    }
    refCnt_++;
    return refCnt_;
}

RefCountIndex Node::decRef() {
    if (refCnt_==0) {
        throw std::length_error("Reference count cannot be negative.");
    }
    refCnt_--;
    return refCnt_;
}

}