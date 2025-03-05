#ifndef __DICTNODEPROVIDER_H__
#define __DICTNODEPROVIDER_H__

#include "INodeProvider.h"
#include "DBManager.h"

class DictNodeProvider : public INodeProvider {
public:
    std::shared_ptr<Node> provide_db_as_node(const std::string& file_path) override;
};

#endif // __DICTNODEPROVIDER_H__
