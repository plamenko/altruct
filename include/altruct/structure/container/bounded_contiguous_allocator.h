#pragma once

#include <memory>

#if _MSC_VER < 1900
#define NOEXCEPT
#else
#define NOEXCEPT noexcept
#endif

namespace altruct {
namespace container {

/**
 * Allocator that provides elements preallocated in a contiguous range.
 *
 * Important: allocator is stateful and hence requires C++11
 * Improtant: only 1 element can be allocated at time
 *
 * @param capacity - maximal number of elements this allocator can hold at any given time
 */
template <class T>
struct bounded_contiguous_allocator {
    typedef T value_type;

    size_t capacity_;
    T* data_;
    size_t available_;
    T** ptrs_;

    bounded_contiguous_allocator(size_t capacity = 1) NOEXCEPT
        : capacity_(capacity), data_(nullptr), available_(capacity), ptrs_(nullptr) {}

    ~bounded_contiguous_allocator() NOEXCEPT{
        clear();
    }

    template <class U>
    bounded_contiguous_allocator(const bounded_contiguous_allocator<U>& allocator) NOEXCEPT
        : bounded_contiguous_allocator(allocator.capacity()) {
    }

    template <class U>
    bounded_contiguous_allocator& operator= (const bounded_contiguous_allocator<U>& allocator) NOEXCEPT{
        clear();
        available_ = capacity_ = allocator.capacity();
        return *this;
    }

    void clear() {
        std::free(ptrs_); ptrs_ = nullptr;
        std::free(data_); data_ = nullptr;
    }

    bool ensure() {
        if (ptrs_) return true;
        if (!data_) data_ = static_cast<T*>(std::malloc(capacity_ * sizeof(T)));
        if (!data_) return false;
        if (!ptrs_) ptrs_ = static_cast<T**>(std::malloc(capacity_ * sizeof(T*)));
        if (!ptrs_) return false;
        for (size_t i = 0; i < capacity_; i++) ptrs_[i] = data_ + i;
        return true;
    }

    size_t capacity() const { return capacity_; }

    T* allocate(size_t n) {
        if (!ensure() || n != 1 || available_ < 1) throw std::bad_alloc();
        return ptrs_[--available_];
    }

    void deallocate(T* ptr, size_t) {
        ptrs_[available_++] = ptr;
    }

    void construct(T* ptr, const T& val) {
        ::new(ptr)T(val);
    }

    void destroy(T* ptr) {
        ptr->~T();
    }
};

} // container
} // altruct

template <class T, class U>
bool operator== (const altruct::container::bounded_contiguous_allocator<T>& lhs,
                 const altruct::container::bounded_contiguous_allocator<U>& rhs) {
    return &lhs == &rhs;
}

template <class T, class U>
bool operator!= (const altruct::container::bounded_contiguous_allocator<T>& lhs,
                 const altruct::container::bounded_contiguous_allocator<U>& rhs) {
    return &lhs != &rhs;
}
