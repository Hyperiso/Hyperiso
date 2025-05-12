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

    if (input_extension == ".slha" || input_extension == ".lha") {
        lha_prototypes.insert(SLHA_BLOCKS.begin(), SLHA_BLOCKS.end());
    } else if (input_extension == ".flha") {
        lha_prototypes.insert(FLHA_BLOCKS.begin(), FLHA_BLOCKS.end());
    }
}   

std::shared_ptr<Node> DBManager::read_from_file(fs::path file_path) {
    // sanitize_file(file_path);  // TODO basic checks (empty file, wrong encoding...)
    ParserFactory::Type parser_type = deduce_parser_type(file_path);
    auto parser = ParserFactory::createParser(parser_type);
    if (parser_type == ParserFactory::Type::LHA) {
        add_default_lha_prototypes(file_path);
        std::dynamic_pointer_cast<LhaParser>(parser)->set_prototypes(this->lha_prototypes);
    }
    auto root = parser->readFromFile(file_path);
    // sanitize_tree(root); // TODO
    return root;
}

void DBManager::write_to_file(fs::path file_path, std::shared_ptr<Node> root) {
    // sanitize_tree(root); // TODO
    ParserFactory::Type parser_type = deduce_parser_type(file_path);
    auto parser = ParserFactory::createParser(parser_type);
    parser->writeToFile(file_path, root);
    // sanitize_file(file_path);  // TODO: basic checks (empty file, wrong encoding...)
}

void DBManager::add_lha_prototype(std::string blockName, size_t itemCount, size_t valueIdx, int scaleIdx, int rgIdx, bool globalScale) {
    std::transform(blockName.begin(), blockName.end(), blockName.begin(), ::toupper);  // Make sure block name is uppercase 
    Prototype new_prototype = Prototype{blockName, itemCount, valueIdx, scaleIdx, rgIdx, globalScale};

    auto it = std::find_if(DBManager::lha_prototypes.begin(), 
                           DBManager::lha_prototypes.end(), 
                           [blockName] (const Prototype& p) { return p.blockName == blockName; });

    if (it != DBManager::lha_prototypes.end())  {
        if (*it == new_prototype) {
            LOG_WARN("Trying to add an already existing prototype for block", blockName);
            return;
        } else {
            LOG_ERROR("ValueError", "Cannot add different prototype for existing block", blockName);
        }
    }

    DBManager::lha_prototypes.emplace(new_prototype);
}
