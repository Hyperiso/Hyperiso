#include "FileWriter.h"

void FileWriter::write(const std::string& dest, std::shared_ptr<BlockAccessor> src) {
    if (src == nullptr) {
        src = MemoryManager::GetInstance()->extract_block_accessor();
    }
    std::shared_ptr<DBNode> node = std::make_shared<DBNode>();

    ParamBlockWriter().write(node, src);

    if (dest.size() >= 5 && dest.substr(dest.size() - 5) == ".json") {
        JSONParser().writeToFile(dest, node);
    } else if (dest.size() >= 5 && dest.substr(dest.size() - 5) == ".yaml") {
        YAMLParser().writeToFile(dest, node);
    } else if (dest.size() >= 5 && dest.substr(dest.size() - 3) == "lha") {
        LhaParser().writeToFile(dest, node);
    } else {

    }
}