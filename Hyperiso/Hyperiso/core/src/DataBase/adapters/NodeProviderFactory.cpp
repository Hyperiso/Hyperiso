#include "NodeProviderFactory.h"

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