#ifndef __NODEPROVIDERFACTORY_H__
#define __NODEPROVIDERFACTORY_H__

#include <memory>
#include "INodeProvider.h"
#include "DictNodeProvider.h"
#include "LhaNodeProvider.h"
#include "DBManager.h"

/**
 * @brief Factory for creating INodeProvider instances.
 */
class NodeProviderFactory {
public:
    enum class Type { DICT, LHA };

    static std::shared_ptr<INodeProvider> createNodeProvider(fs::path src_path);
};

#endif // __NODEPROVIDERFACTORY_H__
