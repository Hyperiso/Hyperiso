#include "DictNodeProvider.h"

std::shared_ptr<Node> DictNodeProvider::provide_db_as_node(const std::string &file_path) {
    return DBManager().read_from_file(file_path);
}