#ifndef __ANSUPPORT_DEFINED
#define __ANSUPPORT_DEFINED

#include <stdint.h>
#include <vector>
#include "common.h"
#include <stdexcept>
#include <iostream>
#include <type_traits>

namespace NAMESPACE {

// We are not supposed to derive from STL classes, so let us 
// rename the template and define some functions. 
template<typename T> using Vector = std::vector<T>;

template<typename T> void zero(std::vector<T>& vec) { 
    vec.assign(vec.size(), 0); 
}

template<typename T> T* dataWithoutBucket(std::vector<T>& vec) { 
    return vec.data()+1; 
}

template<typename T> const T* dataWithoutBucket(const std::vector<T>& vec) { 
    return vec.data()+1; 
}


template<typename T> class CircularBuffer {
public:
    using DepthIndex = HistoryDepthIndex;
    using DepthIndexDelta = HistoryDepthIndexDelta;

    CircularBuffer(DepthIndex size=0) : at_(0), valueCount_(0), size_(size) {
        buffer.resize(size);
    };

    CircularBuffer           (const CircularBuffer&)  = delete;
    CircularBuffer           (      CircularBuffer&&) = delete;
    CircularBuffer& operator=(const CircularBuffer&)  = delete;
    CircularBuffer& operator=(      CircularBuffer&&) = delete;

    // Repository size
    DepthIndex size() const { return size_; };

    // Increase size if new size is greater than current size, otherwise do nothing
    void upsize(DepthIndex newSize) {
        if (newSize<=size_) {
            return;
        }
        buffer.resize(newSize);
        if (size_>0) {
            auto nToMove = size_-at_-1;
            for(auto i=0; i<nToMove; i++) {
                buffer[newSize-1-i] = std::move(buffer[size_-1-i]);
            }
        }
        size_ = newSize;
    };

    // Return the number of values in circular buffer
    DepthIndex valueCount() const { return valueCount_; };

    // Add a value to the buffer
    void add(T value) {
        at_++;
        if (at_>=size_) {
            at_ = 0;
        }
        if (valueCount_<size_) {
            valueCount_++;
        }
        buffer[at_] = value;
    };

    // Remove a value from the buffer
    void remove() {
        if (valueCount_<1) {
            throw std::length_error("Circular buffer is empty. Cannot remove entry.");
        }
        if (at_==0) {
            at_ = size_;
        } else {
            at_--;
        }
        valueCount_--;
    };

    // Get a value from the buffer, throws if which is greater than buffer capacity
    // It is user's responsibility to split the buffer between past and the future
    // 0 = latest value, 1 = second latest value, ... 
    // -1 = first future value, ...
    const T& at(DepthIndexDelta which=0) const {
        if (which>=size_ || -which>=size_) {
            throw std::length_error("Circular buffer is not long enough.");
        }
        if (which>=0) {
            // Past
            if (which>at_) {
                return buffer[size_ + at_ - which];
            } else {
                return buffer[at_ - which];
            }
        } else {
            // Future
            if (at_-which>=size_) {
                return buffer[at_-which-size_]; 
            } else {
                return buffer[at_-which]; 
            }
        }
    };
    
    T& at(DepthIndexDelta which=0) {
        if (which>=size_ || -which>=size_) {
            throw std::length_error("Circular buffer is not long enough.");
        }
        if (which>=0) {
            // Past
            if (which>at_) {
                return buffer[size_ + at_ - which];
            } else {
                return buffer[at_ - which];
            }
        } else {
            // Future
            if (at_-which>=size_) {
                return buffer[at_-which-size_]; 
            } else {
                return buffer[at_-which]; 
            }
        }
    };

    // Advance at by which (>0 = toward future)
    void advance(DepthIndexDelta which=1) {
        if (which>=size_ || -which>=size_) {
            throw std::length_error("Circular buffer is not long enough.");
        }
        if (which<=0) {
            // Past
            if (-which>at_) {
                at_ = size_ + at_ + which;
            } else {
                at_ = at_ + which;
            }
        } else {
            // Future
            if (at_+which>=size_) {
                at_ = at_ + which - size_;
            } else {
                at_ = at_ + which; 
            }
        }
    };

    // Swap at with which (by default swap current and future)
    void rotate(DepthIndexDelta which=-1) {
        if (which>=size_ || -which>=size_) {
            throw std::length_error("Circular buffer is not long enough.");
        }

        DepthIndex atWhich;
        if (which>=0) {
            // Past
            if (which>at_) {
                atWhich = size_ + at_ - which;
            } else {
                atWhich = at_ - which;
            }
        } else {
            // Future
            if (at_-which>=size_) {
                atWhich = at_-which-size_;
            } else {
                atWhich = at_-which; 
            }
        }

        T tmp = std::move(buffer[atWhich]);
        buffer[atWhich] = std::move(buffer[at_]);
        buffer[at_] = std::move(tmp);
    };

    DepthIndex position() const { return at_; }; 

protected:
    DepthIndexDelta size_; // Signed so comparisons with signed will work right
    DepthIndex valueCount_;
    DepthIndex at_;
    std::vector<T> buffer;
};

template<typename T> class VectorRepository : public CircularBuffer<Vector<T>> {
public:
    using DepthIndex = CircularBuffer<std::vector<T>>::DepthIndex;
    using DepthIndexDelta = CircularBuffer<std::vector<T>>::DepthIndexDelta;

    VectorRepository(DepthIndex size=0, size_t length=0);

    VectorRepository           (const VectorRepository&)  = delete;
    VectorRepository           (      VectorRepository&&) = delete;
    VectorRepository& operator=(const VectorRepository&)  = delete;
    VectorRepository& operator=(      VectorRepository&&) = delete;

    // Vector length
    size_t length() const { return length_; };

    // Increase repository size if newSize is greater than size, otherwise do nothing
    void upsize(DepthIndex newSize);

    // Increase repository in terms of size (never decrease), resize all vectors
    void upsize(DepthIndex newSize, size_t newLength);

    // Set a vector to 0.0, which=0 corresponds to current vector, 
    // which>0 corresponds to past vectors
    // which<0 corresponds to future vectors
    void zero(DepthIndex which=0);
    void zeroPast() { zero(1); };
    void zeroFuture() { zero(-1); };

    using CircularBuffer<Vector<T>>::advance;

    using CircularBuffer<Vector<T>>::at;
    
    using CircularBuffer<Vector<T>>::rotate;

    // which=0 returns current vector
    // which>0 returns past vector
    // which<0 returns future vector
    // abs(which)>=size throws
    
    // Returns vector of doubles including bucket
    Vector<T>& vector(DepthIndex which=0) { return at(which); };
    const Vector<T>& vector(DepthIndex which=0) const { return at(which); };
    Vector<T>& pastVector() { return at(1); };
    const Vector<T>& pastVector() const { return at(1); };
    Vector<T>& futureVector() { return at(-1); };
    const Vector<T>& futureVector() const { return at(-1); };
    
    // Return array of doubles including bucket
    T* data(DepthIndexDelta which=0) { return at(which).data(); };
    const T* data(DepthIndexDelta which=0) const { return at(which).data(); };
    T* pastData() { return data(1); };
    const T* pastData() const { return data(1); };
    T* futureData() { return data(-1); };
    const T* futureData() const { return data(-1); };
    
    // Return array of doubles without bucket 
    // (for use with linear solver)
    T* dataWithoutBucket(DepthIndexDelta which=0) { return dataWithoutBucket(vector(which)); };
    const T* dataWithoutBucket(DepthIndexDelta which=0) const { return dataWithoutBucket(vector(which)); }
    T* pastDataWithoutBucket() { return dataWithoutBucket(1); };
    const T* pastDataWithoutBucket() const { return dataWithoutBucket(1); }
    T* futureDataWithoutBucket() { return dataWithoutBucket(-1); };
    const T* futureDataWithoutBucket() const { return dataWithoutBucket(-1); }

protected:
    using CircularBuffer<Vector<T>>::size_;
    using CircularBuffer<Vector<T>>::valueCount_;
    using CircularBuffer<Vector<T>>::at_;
    using CircularBuffer<Vector<T>>::buffer;

private:
    size_t length_;
};

// Template implementation (verbose members)
template<typename T> VectorRepository<T>::VectorRepository(DepthIndex size, size_t length) 
    : CircularBuffer<Vector<T>>(size), length_(length) {
    for(auto i=0; i<size_; i++) {
        buffer[i].resize(length_); 
    }
}

template<typename T> void VectorRepository<T>::upsize(DepthIndex newSize) {
    if (newSize<=size_) {
        return;
    }
    auto oldSize = size_;
    CircularBuffer<Vector<T>>::upsize(newSize);
    // Extend new vectors
    auto delta = newSize - oldSize;
    for(decltype(delta) i=0; i<delta; i++) {
        CircularBuffer<Vector<T>>::at(newSize-1-i).resize(length_);
    }
}

template<typename T> void VectorRepository<T>::upsize(DepthIndex newSize, size_t newLength) {
    for(auto i=0; i<size_; i++) {
        buffer[i].resize(newLength); 
    }
    length_ = newLength;
    upsize(newSize);
    return;
}

template<typename T> void VectorRepository<T>::zero(DepthIndex which) {
    if (length_>0) {
        CircularBuffer<Vector<T>>::at(which).assign(length_, 0.0);
    }
}

}

#endif
