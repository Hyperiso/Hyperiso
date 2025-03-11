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
    
    std::shared_ptr<INodeProvider> NodeProviderFactory::createNodeProvider(fs::path src_path) {
        auto input_extension = src_path.extension().string();
        auto dict_extensions = DBManager::EXTENSIONS.at(ParserFactory::Type::JSON);
        dict_extensions.insert(DBManager::EXTENSIONS.at(ParserFactory::Type::YAML).begin(), DBManager::EXTENSIONS.at(ParserFactory::Type::YAML).end());
        auto lha_extensions = DBManager::EXTENSIONS.at(ParserFactory::Type::LHA);

        if (dict_extensions.contains(input_extension)) {
            return std::make_shared<DictNodeProvider>(src_path);
        } else if (lha_extensions.contains(input_extension)) {
            return std::make_shared<LhaNodeProvider>(src_path);
        }

        LOG_ERROR("FileError", "Cannot read file with extension", input_extension);
    }

#endif // __NODEPROVIDERFACTORY_H__
