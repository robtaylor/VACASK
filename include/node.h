#ifndef __NODE_DEFINED
#define __NODE_DEFINED

#include <type_traits>
#include "flags.h"
#include "identifier.h"
#include "common.h"


namespace NAMESPACE {

enum class NodeFlags : uint8_t {
    NodeTypeMask = 1,
    PotentialNode = 0,
    FlowNode = 1,

    InternalDeviceNode = 2,
    Ground = 4,
};
DEFINE_FLAG_OPERATORS(NodeFlags);

class Node : public FlagBase<NodeFlags> {
public:
    Node(Id name, Flags f);
    ~Node();

    Node           (const Node&)  = delete;
    Node           (      Node&&) = delete;
    Node& operator=(const Node&)  = delete;
    Node& operator=(      Node&&) = delete;

    Id name() { return name_; };
    
    void setUnknownIndex(UnknownIndex u) { unknown_=u; };
    UnknownIndex unknownIndex() { return unknown_; };

    RefCountIndex incRef();
    RefCountIndex decRef();
    
private:
    Id name_;
    Flags flags_;
    UnknownIndex unknown_;
    RefCountIndex refCnt_;
};

}

#endif
