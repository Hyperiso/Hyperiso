#include "DBManager.h"

ParserFactory::Type DBManager::deduce_parser_type(fs::path file_path) {
    auto input_extension = file_path.extension().string();
    for (const auto& allowed_exts : DBManager::EXTENSIONS) {
        if (allowed_exts.second.contains(input_extension))
            return allowed_exts.first;
    }

    LOG_ERROR("FileError", "Cannot read file with extension", input_extension);
}

void DBManager::add_default_lha_prototypes(fs::path file_path) {
    auto input_extension = file_path.extension().string();
    lha_prototypes.insert(LHA_BLOCKS.begin(), LHA_BLOCKS.end());

    if (input_extension == ".slha") {
        lha_prototypes.insert(SLHA_BLOCKS.begin(), SLHA_BLOCKS.end());
    } else if (input_extension == ".flha") {
        lha_prototypes.insert(FLHA_BLOCKS.begin(), FLHA_BLOCKS.end());
    }
}   

std::shared_ptr<Node> DBManager::read_from_file(fs::path file_path) {
    ParserFactory::Type parser_type = deduce_parser_type(file_path);
    // sanitize_file(file_path);  // basic checks (empty file, wrong encoding...)
    auto parser = ParserFactory::createParser(parser_type);
    if (parser_type == ParserFactory::Type::LHA) {
        add_default_lha_prototypes(file_path);
        static_cast<std::shared_ptr<LhaParser>>(parser)->setPrototypes(lha_prototypes);
    }
    auto root = parser->readFromFile(file_path);
    // sanitize_tree(root);
    return root;
}

void DBManager::write_to_file(fs::path file_path, std::shared_ptr<Node> root) {
    // sanitize_tree(root);
    ParserFactory::Type parser_type = deduce_parser_type(file_path);
    auto parser = ParserFactory::createParser(parser_type);
    parser->writeToFile(file_path, root);
    // sanitize_file(file_path);  // basic checks (empty file, wrong encoding...)
}

void DBManager::add_lha_prototype(std::string blockName, int itemCount=2, int valueIdx=1, int scaleIdx=-1, int rgIdx=-1, bool globalScale=false) {
    std::transform(blockName.begin(), blockName.end(), blockName.begin(), ::toupper);  // Make sure block name is uppercase 
    this->lha_prototypes.emplace(Prototype{blockName, itemCount, valueIdx, scaleIdx, rgIdx, globalScale});
}
