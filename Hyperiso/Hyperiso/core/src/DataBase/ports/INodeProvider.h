#ifndef __INODEPROVIDER_H__
#define __INODEPROVIDER_H__

#include <memory>
#include "DBNode.h"

class INodeProvider {
public:
    virtual ~INodeProvider() = default;

    virtual std::shared_ptr<Node> provide_db_as_node(const std::string& file_path) = 0; 
};

#endif // __INODEPROVIDER_H__
