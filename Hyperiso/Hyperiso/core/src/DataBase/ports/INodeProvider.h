#ifndef __INODEPROVIDER_H__
#define __INODEPROVIDER_H__

#include <memory>
#include <filesystem>
#include "DBNode.h"

namespace fs = std::filesystem;

class INodeProvider {
protected:
    fs::path src_path;

public:
    INodeProvider(fs::path src_path) : src_path(src_path) {}
    virtual ~INodeProvider() = default;

    virtual std::shared_ptr<Node> provide_db_as_node() = 0; 
};

#endif // __INODEPROVIDER_H__
