#include <iostream>
#include <unordered_map>
#include <memory>
#include <cstdlib>

class MemoryManager {
private:
    DataBase* db;
    std::unordered_map<int, std::unique_ptr<std::string>> cache;

public:
    explicit MemoryManager(DataBase* db) : db(db) {}

    // Creation pointeur unique
    template<typename T>
    std::unique_ptr<T, void(*)(void*)> makeUniquePtr(T* ptr) {
        return std::unique_ptr<T, void(*)(void*)>(ptr, [](void* p) { std::free(p); });
    }

    // Méthode pour allouer de la mémoire
    template<typename T>
    std::unique_ptr<T, void(*)(void*)> allocate() {
        T* ptr = static_cast<T*>(std::malloc(sizeof(T)));
        if (!ptr) {
            throw std::bad_alloc();
        }
        return makeUniquePtr(ptr);
    }

    const std::string& getData(int key) {
        if (cache.find(key) == cache.end()) {
            cache[key] = std::make_unique<std::string>(db->getData(key));
        }
        return *cache[key];
    }

    // Autres méthodes...

    MemoryManager(const MemoryManager&) = delete;
    MemoryManager& operator=(const MemoryManager&) = delete;
    MemoryManager(MemoryManager&&) noexcept = default;
    MemoryManager& operator=(MemoryManager&&) noexcept = default;
};

