#ifndef __BLOCKSCREATOR_H__
#define __BLOCKSCREATOR_H__

#include <memory>
#include "BlockAccessor.h"
#include "lha_reader.h"
#include "DBNode.h"

class BlocksCreator {
public:
    static std::shared_ptr<BlockAccessor> from_lha_reader(std::shared_ptr<LhaReader> reader);
    static std::shared_ptr<BlockAccessor> from_db_node(std::shared_ptr<Node> root);
};

#endif // __BLOCKSCREATOR_H__
