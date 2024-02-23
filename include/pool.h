#ifndef __POOL_DEFINED
#define __POOL_DEFINED

#include <forward_list>
#include <iostream>
#include <cstring>
#include <vector>
#include <utility>
#include "common.h"


namespace NAMESPACE {

// Pool allocator, returns pointer of type void* to allocated space
class PoolAllocator {
public:
    typedef struct Block {
        void* mem;
        size_t used;
        size_t size;
        size_t failures;

        // Disable default constructor
        Block() = delete;

        Block(size_t size) : size(size), used(0), failures(0) {
            mem = malloc(size);
        }; 

        ~Block() {
            free(mem);
        }

        void* basePtr(size_t offset=0) { return reinterpret_cast<void*>(reinterpret_cast<char*>(mem)+offset); };
        
        void* allocate(size_t count) {
            if (used+count<=size) {
                auto ptr = basePtr(used);
                used += count;
                return ptr;
            }
            return nullptr;
        }
    } Block;

    // Disable default constructor
    PoolAllocator() = delete;

    PoolAllocator(size_t entrySize, size_t initialCount=1024, size_t growthFactor=2, size_t maxFailures=4);
    virtual ~PoolAllocator();
    virtual void clear();
    void* allocate(size_t count);
    void dump(std::ostream& os);

protected:
    using BlockList = std::forward_list<Block>;
    BlockList blockList;
    void advance();
    using BlockIterator = typename BlockList::iterator;
    BlockIterator atBlock;
    BlockIterator lastBlock;
    size_t initialCount;
    size_t entrySize;
    size_t growthFactor;
    size_t maxFailures;
};

// Pool allocator with entries of type T, returns pointer of type T* to allocated space
template<typename T, bool doPlacement=true> class TypedPoolAllocator : public PoolAllocator {
public:
    TypedPoolAllocator(size_t initialCount=1024, size_t growthFactor=2, size_t maxFailures=4);
    virtual ~TypedPoolAllocator();
    virtual void clear();
    template<class... Args> T* allocate(Args&&... args);
    template<class... Args> T* allocate(size_t count, Args&&... args);
};

template<typename T, bool doPlacement> TypedPoolAllocator<T, doPlacement>::TypedPoolAllocator(size_t initialCount, size_t growthFactor, size_t maxFailures) 
    : PoolAllocator(sizeof(T), initialCount, growthFactor, maxFailures) {
}

template<typename T, bool doPlacement> TypedPoolAllocator<T, doPlacement>::~TypedPoolAllocator() {
    // Call destructor on entries if we used placement
    if constexpr(doPlacement) {
        for(auto& block : blockList) {
            auto end = reinterpret_cast<T*>(block.basePtr(block.used));
            for(auto ptr = reinterpret_cast<T*>(block.basePtr()); ptr<end; ++ptr) {
                ptr->~T();
            }
        }
    }
}

template<typename T, bool doPlacement> void TypedPoolAllocator<T, doPlacement>::clear() {
    // Call destructor on entries if we used placement
    if constexpr(doPlacement) {
        for(auto& block : blockList) {
            auto end = reinterpret_cast<T*>(block.basePtr(block.used));
            for(auto ptr = reinterpret_cast<T*>(block.basePtr()); ptr<end; ++ptr) {
                ptr->~T();
            }
        }
    }
    // Free blocks
    PoolAllocator::clear();
}

template<typename T, bool doPlacement> template<class... Args> T* TypedPoolAllocator<T, doPlacement>::allocate(Args&&... args) {
    auto ptr = PoolAllocator::allocate(1);
    auto objPtr = reinterpret_cast<T*>(ptr);
    // Placement new
    if constexpr(doPlacement) {
        // Placement new with custom constructor
        return new (objPtr) T(std::forward<Args>(args)...);
    } else {
        // No placement new
        return objPtr;
    }
}

template<typename T, bool doPlacement> template<class... Args> T* TypedPoolAllocator<T, doPlacement>::allocate(size_t count, Args&&... args) {
    auto ptr = PoolAllocator::allocate(count);
    auto objPtr = reinterpret_cast<T*>(ptr);
    if constexpr(doPlacement) {
        for(decltype(count) i=0; i<count; i++) {
            // Placement new with default constructor
            new (objPtr+i) T();
        }
    } 
    return objPtr;
}


// Pool allocator with entries of type T, entries can be freed
// Upon new allocation freed entries are recycled
// allocate() returns only a single entry
template<typename T, bool doPlacement=true> class RecyclingTypedPoolAllocator : private TypedPoolAllocator<T, doPlacement> {
public:
    RecyclingTypedPoolAllocator(size_t initialCount=1024, size_t growthFactor=2, size_t maxFailures=4);
    virtual ~RecyclingTypedPoolAllocator();
    virtual void clear();
    
    template<class... Args> T* allocate(Args&&... args);
    void free(T* ptr);

private:
    std::vector<T*> freeEntries;
};

template<typename T, bool doPlacement> RecyclingTypedPoolAllocator<T, doPlacement>::RecyclingTypedPoolAllocator(size_t initialCount, size_t growthFactor, size_t maxFailures) 
    : TypedPoolAllocator<T, doPlacement>(initialCount, growthFactor, maxFailures) {
}

template<typename T, bool doPlacement> RecyclingTypedPoolAllocator<T, doPlacement>::~RecyclingTypedPoolAllocator() {
}

template<typename T, bool doPlacement> void RecyclingTypedPoolAllocator<T, doPlacement>::clear() {
    // Call destructors if placement was used, free blocks
    TypedPoolAllocator<T, doPlacement>::clear();
    // Clear the list of freed objects
    freeEntries.clear();
}

template<typename T, bool doPlacement> template<class... Args> T* RecyclingTypedPoolAllocator<T, doPlacement>::allocate(Args&&... args) {
    if (freeEntries.size()>0) {
        auto ptr = freeEntries.back();
        freeEntries.pop_back();
        if constexpr(doPlacement) {
            ptr->~T();
            new (ptr) T(std::forward<Args>(args)...);
        }
        return ptr;
    }
    return TypedPoolAllocator<T>::allocate(std::forward<Args>(args)...);
}

template<typename T, bool doPlacement> void RecyclingTypedPoolAllocator<T, doPlacement>::free(T* ptr) {
    // Do not call the destructor, do it when the entry is recycled
    freeEntries.push_back(ptr);
}


// Pool of strings, strings cannot be freed individually
class CStringPool : public TypedPoolAllocator<char, false> {
public:
    CStringPool(size_t blockSize=8192, size_t growthFactor=2, size_t maxFailures=4);
    virtual ~CStringPool();

    const char* allocate(const char* s);
    const char* allocate(const char* s, size_t n);
    const char* allocate(const std::string& s);

    void dump(std::ostream& os, bool details=true);
};

}

#endif
