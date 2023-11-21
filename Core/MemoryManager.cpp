#include <cstddef>
#include <iostream>
#include <unordered_map>

class MemoryManager {
public:
    static MemoryManager& getInstance() {
        static MemoryManager instance;
        return instance;
    }

    void* allocate(std::size_t size) {
        void* ptr = std::malloc(size);
        allocations[ptr] = size;
        std::cout << "Allocated " << size << " bytes at " << ptr << std::endl;
        return ptr;
    }

    void deallocate(void* ptr) {
        if (allocations.find(ptr) != allocations.end()) {
            std::cout << "Deallocated " << allocations[ptr] << " bytes from " << ptr << std::endl;
            std::free(ptr);
            allocations.erase(ptr);
        } else {
            std::cout << "Attempt to deallocate unrecognized memory: " << ptr << std::endl;
        }
    }

    // Add more functionality as needed

private:
    MemoryManager() {}
    ~MemoryManager() {
        // Cleanup code if required
    }

    std::unordered_map<void*, std::size_t> allocations;
};
