#ifndef __GENERATOR_DEFINED
#define __GENERATOR_DEFINED

#include <concepts>
#include <coroutine>
#include <exception>
#include "common.h"

// Based on: https://www.scs.stanford.edu/~dm/blog/c++-coroutines.html

namespace NAMESPACE {

template<typename T>
struct Generator {
    // Promise type and handle
    struct promise_type;
    using handle_type = std::coroutine_handle<promise_type>;

    // Promise structure
    struct promise_type {
        // Value storage
        T value_;

        // Exception that happened in the coroutine
        std::exception_ptr exception_;

        // Create generator
        // When a coroutine begins execution, the return value is stored and 
        // returned to the caller when the coroutine suspends for the first time. 
        Generator get_return_object() {
            return Generator(handle_type::from_promise(*this));
        }

        // At the beginning of the coroutine there is an implicit co_await
        // The return value tells what it should do. 
        // suspend_always means that the coroutine will produce 
        // the first value at the first suspension in its body
        std::suspend_always initial_suspend() { 
            return {}; 
        }
        
        // At the end of the coroutine is an implicit co_await
        // The return value tells what it should do
        std::suspend_always final_suspend() noexcept { 
            return {}; 
        }
        
        // Called from implicit catch block when an unhandled exception in coroutine happens
        void unhandled_exception() { 
            exception_ = std::current_exception(); 
        }
        
        // Called when co_yield is invoked
        template<std::convertible_to<T> From> // C++20 concept
        std::suspend_always yield_value(From&& from) {
            value_ = std::forward<From>(from);
            return {};
        }
        
        // Called when co_return with no arguments is invoked
        void return_void() {
        }

        // Called when co_return with an arguments is invoked
        // We cannot have both return_void() and return_value()
        // template<std::convertible_to<T> From> // C++20 concept
        //     void return_value(From&& from) {
        // }
    };

    // Promise handle
    handle_type h_;

    // Default constructor
    Generator() : h_(nullptr) {}; 

    // Copy constructor is disabled
    Generator(const Generator&) = delete;

    // Move constructor
    Generator(Generator&& other) {
        if (h_!=nullptr) {
            h_.destroy();
        }
        h_ = std::exchange(other.h_, nullptr);
    };

    // Copy assignment is disabled
    Generator& operator=(const Generator&)  = delete;
    
    // Move assignment
    Generator& operator=(Generator&& other) {
        if (h_!=nullptr) {
            h_.destroy();
        }
        h_ = std::exchange(other.h_, nullptr);
        return *this;
    }
    
    // Cleanup - destroy coroutine
    ~Generator() { 
        if (h_!=nullptr) {
            h_.destroy(); 
            h_ = nullptr;
        }
    }

    bool isValid() const { return h_!=nullptr; };

    // Check if the coroutine is not finished yet, must be a valid coroutine
    bool done() {
        // Fill with new value if needed
        fill();
        // Check if the coroutine is done
        // done() returns true if the coroutine is suspended at its final suspend point
        return h_.done();
    } 

    // Convert generator to bool for checking if the coroutine is not finished yet
    // Same as done()
    explicit operator bool() {
        return !done();
    }
    
    // Generator call, proceed with coroutine if no value available, return yielded value
    T resume() {
        // Fill promise with value, if one is not already there
        fill();
        // Remove the value and return it
        full_ = false;
        return std::move(h_.promise().value_);
    }

    // Same as operator()
    T operator()() {
        return std::move(resume());
    };
    
private:
    // Constructor
    Generator(handle_type h) : h_(h) {}

    // Do we have a value in the promise? 
    bool full_ = false;

    // Fill promise with value
    void fill() {
        // Value already there? 
        if (!full_) {
            // No, resume coroutine
            h_.resume();
            // In case of exception rethrow
            if (h_.promise().exception_)
                std::rethrow_exception(h_.promise().exception_);
            // We have a value now
            full_ = true;
        }
    };
};

/*
// Example

Generator<int> counter(int from, int to) {
    for(int i=from; i<to; i++) {
        co_yield i;
    }
}

int main() {
    // Create a generator, run it until initial suspend
    auto gen = counter(0, 10);

    while (gen) {
        std::cout << gen() << "\n";
    }

    return 0;
}
*/

}

#endif
