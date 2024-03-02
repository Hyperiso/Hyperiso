#include "MemoryManager.h"
#include <iostream>
#include <unordered_map>
#include <memory>
#include <cstdlib>
#include "../DataBase/lha_reader.h"

// Creation pointeur unique
template<typename T>
std::unique_ptr<T, void(*)(void*)> MemoryManager::makeUniquePtr(T* ptr) {
    return std::unique_ptr<T, void(*)(void*)>(ptr, [](void* p) { std::free(p); });
}

// Méthode pour allouer de la mémoire
template<typename T>
std::unique_ptr<T, void(*)(void*)> MemoryManager::allocate() {
    T* ptr = static_cast<T*>(std::malloc(sizeof(T)));
    if (!ptr) {
        throw std::bad_alloc();
    }
    return makeUniquePtr(ptr);
}

MemoryManager* MemoryManager::GetInstance(std::string lhaFile="./example.flha", std::vector<int> models={0}) {
    if (!instance) {
        MemoryManager::instance = new MemoryManager(lhaFile, models);
    }
    return MemoryManager::instance;
}

const std::string &MemoryManager::getData(std::string key)
{
    // if (cache.find(key) == cache.end()) {
    //     cache[key] = std::make_unique<std::string>(db->readBlock(key));
    // }
    return *cache[key];
}

LhaReader* MemoryManager::getReader() {
    return reader.get();
}

void MemoryManager::init() {
    // Init SM parameters from SMINPUTS
    reader = std::make_unique<LhaReader>(LhaReader(lhaPath));
    reader->readAll();

    // Init other parameters

}
