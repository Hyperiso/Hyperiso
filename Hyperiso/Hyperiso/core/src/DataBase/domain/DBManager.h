#ifndef __DBMANAGER_H__
#define __DBMANAGER_H__

#include <memory>
#include <filesystem>
#include "DBNode.h"
#include "ParserFactory.h"
#include "LhaBlockPrototype.h"

namespace fs = std::filesystem;

class DBManager {
private:
    static inline const std::map<ParserFactory::Type, std::unordered_set<std::string>> EXTENSIONS {
        {ParserFactory::Type::JSON, {".json"}},
        {ParserFactory::Type::YAML, {".yml", ".yaml"}},
        {ParserFactory::Type::LHA, {".lha", ".slha", ".flha"}},
    }; 

    std::unordered_set<Prototype> lha_prototypes;

    ParserFactory::Type deduce_parser_type(fs::path file_path);

    void add_default_lha_prototypes(fs::path file_path);

public:
    std::shared_ptr<Node> read_from_file(fs::path file_path);
    void write_to_file(fs::path file_path, std::shared_ptr<Node> root);

    /**
     * @brief Adds a new block type prototype to the reader.
     * @param blockName Name of the block.
     * @param itemCount Number of items in the block.
     * @param valueIdx Index of the value column.
     * @param scaleIdx Index of the scale column.
     * @param rgIdx Index of the renormalization group column.
     * @param globalScale Flag indicating if the block uses a global scale.
     */
    void add_lha_prototype(std::string blockName, int itemCount=2, int valueIdx=1, int scaleIdx=-1, int rgIdx=-1, bool globalScale=false);
};


#endif // __DBMANAGER_H__
