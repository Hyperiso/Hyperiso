#ifndef __DICTNODEPROVIDER_H__
#define __DICTNODEPROVIDER_H__

#include "INodeProvider.h"
#include "DBManager.h"

class DictNodeProvider : public INodeProvider {
public:
    DictNodeProvider(fs::path src_path) : INodeProvider(src_path) {}

    std::shared_ptr<Node> provide_db_as_node() override;
};

#endif // __DICTNODEPROVIDER_H__
