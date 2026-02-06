#include "DictDBNodeProvider.h"

std::shared_ptr<DBNode> DictDBNodeProvider::provide_db_as_node() {
    return DBManager().read_from_file(this->src_path);
}