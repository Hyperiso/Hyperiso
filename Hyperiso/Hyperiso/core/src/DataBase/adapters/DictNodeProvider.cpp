#include "DictNodeProvider.h"

std::shared_ptr<Node> DictNodeProvider::provide_db_as_node() {
    return DBManager().read_from_file(this->src_path);
}