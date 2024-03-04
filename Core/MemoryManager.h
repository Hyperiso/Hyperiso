#if !defined(HYPERISO_MEMORY_MANAGER_H)
#define HYPERISO_MEMORY_MANAGER_H

#include <iostream>
#include <unordered_map>
#include <memory>
#include <cstdlib>
#include "lha_reader.h"

class MemoryManager {
private:
    std::string lhaPath;
    std::vector<int> models;
    std::unordered_map<std::string, std::unique_ptr<std::string>> cache;
    std::unique_ptr<LhaReader> reader;
    static MemoryManager* instance;
    inline explicit MemoryManager(std::string lhaPath, std::vector<int> models) : lhaPath(lhaPath), models(std::move(models)) {};

public:

    static MemoryManager* GetInstance(std::string lhaFile="./example.flha", std::vector<int> models={0});

    // Creation pointeur unique
    template<typename T>
    std::unique_ptr<T, void(*)(void*)> makeUniquePtr(T* ptr);

    // Méthode pour allouer de la mémoire
    template<typename T>
    std::unique_ptr<T, void(*)(void*)> allocate();

    const std::string& getData(std::string key);

    LhaReader* getReader();

    // initializes the memory with all the necessary parameters and those read in the LHA file
    void init();

    // Autres méthodes...

    MemoryManager(const MemoryManager&) = delete;
    MemoryManager& operator=(const MemoryManager&) = delete;
    MemoryManager(MemoryManager&&) noexcept = default;
    MemoryManager& operator=(MemoryManager&&) noexcept = default;
};



#endif // HYPERISO_MEMORY_MANAGER_H
