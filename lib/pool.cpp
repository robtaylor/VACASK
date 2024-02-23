#include "pool.h"

namespace NAMESPACE {

PoolAllocator::PoolAllocator(size_t entrySize, size_t initialCount, size_t growthFactor, size_t maxFailures) 
    : entrySize(entrySize), initialCount(initialCount), growthFactor(growthFactor), maxFailures(maxFailures) {
    blockList.emplace_front(initialCount*entrySize);
    atBlock = blockList.begin();
    lastBlock = atBlock;
}

PoolAllocator::~PoolAllocator() {
    blockList.clear();
}

void PoolAllocator::clear() {
    blockList.clear();
    blockList.emplace_front(initialCount*entrySize);
    atBlock = blockList.begin();
    lastBlock = atBlock;
}

void* PoolAllocator::allocate(size_t count) {
    // Compute needed size
    auto need = count*entrySize;
    // Scan blocks for free space
    for(auto at=atBlock; at!=blockList.end(); ++at) {
        // Try allocation
        auto ptr = at->allocate(need);
        if (ptr) {
            // Found a block with sufficient free space
            // Skip to first open block pointer
            advance();
            return ptr;
        }
        // Increase failure count for examined block
        at->failures++;
    }
    // Not found, create new block
    // Try last size times growthFactor
    // If not enough, use growthFactor times size of allocation
    size_t newSize = lastBlock->size * growthFactor;
    if (newSize<need) {
        newSize = need * growthFactor;
    }
    auto newBlock = blockList.emplace_after(lastBlock, newSize);
    ++lastBlock;
    auto newPtr = newBlock->allocate(need);
    return newPtr;
}

void PoolAllocator::advance() {
    for(auto at = atBlock; at!=blockList.end(); ++at) {
        // Stop at first block with failures below threshold
        if (at->failures<maxFailures) {
            atBlock = at;
            return;
        }
    }
    // Nothing found, move to last block
    atBlock = lastBlock;
}

void PoolAllocator::dump(std::ostream& os) {
    for(auto at=blockList.begin(); at!=blockList.end(); ++at) {
        os << "Block size=" << at->size << ", used=" << at->used << std::endl; 
    }
}


CStringPool::CStringPool(size_t blockSize, size_t growthFactor, size_t maxFailures) 
    : TypedPoolAllocator<char, false>(blockSize, growthFactor, maxFailures) {
}

CStringPool::~CStringPool() {
}

const char* CStringPool::allocate(const char* s, size_t n) {
    auto ptr = TypedPoolAllocator<char, false>::allocate(n+1);
    std::strcpy(ptr, s);
    return ptr;
}

const char* CStringPool::allocate(const char* s) {
    return allocate(s, std::strlen(s));
}

const char* CStringPool::allocate(const std::string& s) {
    return allocate(s.c_str(), s.size());
}

void CStringPool::dump(std::ostream& os, bool details) {
    for(auto at=blockList.begin(); at!=blockList.end(); ++at) {
        os << "Block size " << at->size << ", used " << at->used << std::endl; 
        size_t strCount = 0;
        for(size_t i=0; i<at->used; ) {
            auto ptr = reinterpret_cast<char*>(at->basePtr(i));
            auto n = std::strlen(ptr);
            if (details) {
                os << "  " << strCount << ":" << ptr << std::endl;
            }
            i += n+1;
            strCount++;
        }
        std::cout << "  string count: " << strCount << std::endl;
    }
}

}
