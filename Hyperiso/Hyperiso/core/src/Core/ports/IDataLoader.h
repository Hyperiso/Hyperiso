#ifndef __IDATALOADER_H__
#define __IDATALOADER_H__

#include <memory>
#include <filesystem>

namespace fs = std::filesystem;

template<typename T>
class IDataLoader {
public:
    virtual void load(std::shared_ptr<T> dest, fs::path src_file) = 0;

    virtual ~IDataLoader() = default;
};

#endif // __IDATALOADER_H__
