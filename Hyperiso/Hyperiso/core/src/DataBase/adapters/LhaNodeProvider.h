#ifndef __LHANODEPROVIDER_H__
#define __LHANODEPROVIDER_H__

#include <string>
#include "INodeProvider.h"
#include "DBManager.h"

class LhaNodeProvider : public INodeProvider {
public:
    LhaNodeProvider(fs::path src_path) : INodeProvider(src_path) {}

    std::shared_ptr<Node> provide_db_as_node() override;

    /**
     * @brief Adds a new block type prototype to the reader.
     * @param blockName Name of the block.
     * @param itemCount Number of items in the block.
     * @param valueIdx Index of the value column.
     * @param scaleIdx Index of the scale column.
     * @param rgIdx Index of the renormalization group column.
     * @param globalScale Flag indicating if the block uses a global scale.
     */
    void add_lha_prototype(BlockName blockName, size_t itemCount=2, size_t valueIdx=1, int scaleIdx=-1, int rgIdx=-1, bool globalScale=false);

};

#endif // __LHANODEPROVIDER_H__
