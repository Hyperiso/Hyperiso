#ifndef __IPATHPROVIDER_H__
#define __IPATHPROVIDER_H__

#include <filesystem>

namespace fs = std::filesystem;

template <typename EnumType>
class IPathProvider {
public:
    virtual ~IPathProvider() = default;

    virtual fs::path get_path(EnumType);
};


#endif // __IPATHPROVIDER_H__
