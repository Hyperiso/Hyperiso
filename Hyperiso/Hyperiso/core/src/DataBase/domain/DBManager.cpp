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

    auto insert_if_absent = [&](const Prototype& p) {
        auto it = std::find_if(lha_prototypes.begin(), lha_prototypes.end(),
            [&](const Prototype& q){ return q.blockName == p.blockName; });
        if (it == lha_prototypes.end()) {
            lha_prototypes.insert(p);
        }
    };

    for (const auto& p : LHA_BLOCKS) insert_if_absent(p);

    if (input_extension == ".slha") {
        for (const auto& p : SLHA_BLOCKS) insert_if_absent(p);
    } else if (input_extension == ".flha") {
        for (const auto& p : FLHA_BLOCKS) insert_if_absent(p);
    } else if (input_extension == ".lha") {
        // The generic .lha extension is used by both SLHA-like spectrum
        // generators (including 2HDMC) and FLHA files.  Accept both
        // standard prototype sets so valid mixed-format metadata such as
        // FMODSEL does not generate unsupported-block warnings.
        for (const auto& p : SLHA_BLOCKS) insert_if_absent(p);
        for (const auto& p : FLHA_BLOCKS) insert_if_absent(p);
    }
} 

void DBManager::sanitize_file(const fs::path& file_path) {
    if (!fs::exists(file_path)) {
        throw std::runtime_error("File does not exist: " + file_path.string());
    }

    if (fs::is_empty(file_path)) {
        throw std::runtime_error("File is empty: " + file_path.string());
    }

    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + file_path.string());
    }

    char c;
    while (file.get(c)) {
        if (static_cast<unsigned char>(c) < 0x09 && c != '\n' && c != '\r' && c != '\t') {
            throw std::runtime_error("File may contain invalid characters.");
        }
    }
}

std::shared_ptr<DBNode> DBManager::read_from_file(fs::path file_path) {
    sanitize_file(file_path);
    ParserFactory::Type parser_type = deduce_parser_type(file_path);
    auto parser = ParserFactory::createParser(parser_type);
    if (parser_type == ParserFactory::Type::LHA) {
        add_default_lha_prototypes(file_path);
        std::dynamic_pointer_cast<LhaParser>(parser)->set_prototypes(this->lha_prototypes);
    }
    auto root = parser->readFromFile(file_path);
    sanitize_tree(root);
    return root;
}

void DBManager::sanitize_tree(const std::shared_ptr<DBNode>& root, const std::vector<BlockName>& required_keys) {
    if (!root) {
        throw std::runtime_error("Tree is null.");
    }

    if (root->countChildren() == 0) {
        throw std::runtime_error("Tree is empty.");
    }

    for (const auto& key : required_keys) {
        if (!root->contains(key)) {
            throw std::runtime_error("Missing required key in tree: " + key);
        }
    }
}


void DBManager::write_to_file(fs::path file_path, std::shared_ptr<DBNode> root) {
    sanitize_tree(root, { "SMINPUTS" });
    ParserFactory::Type parser_type = deduce_parser_type(file_path);
    auto parser = ParserFactory::createParser(parser_type);
    parser->writeToFile(file_path, root);
    sanitize_file(file_path);
}

void DBManager::add_lha_prototype(BlockName blockName, size_t itemCount, size_t valueIdx, int scaleIdx, int rgIdx, int binIdx, bool globalScale) {
    blockName.to_upper();
    Prototype new_prototype = Prototype{blockName, itemCount, valueIdx, scaleIdx, rgIdx, binIdx, globalScale};

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
