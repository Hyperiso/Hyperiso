#include "LhaNodeProvider.h"

std::shared_ptr<Node> LhaNodeProvider::provide_db_as_node() {
    return DBManager().read_from_file(this->src_path);
}

void LhaNodeProvider::add_lha_prototype(BlockName blockName,
                                        size_t itemCount,
                                        size_t valueIdx,
                                        int scaleIdx,
                                        int rgIdx,
                                        int binIdx,
                                        bool globalScale)
{
    bool input_ok = valueIdx < itemCount && scaleIdx < itemCount && rgIdx < itemCount && binIdx + 1 < itemCount
                    && scaleIdx >= -1 && rgIdx >= -1 && binIdx >= -1 
                    && valueIdx != scaleIdx && valueIdx != rgIdx && valueIdx != binIdx && valueIdx != binIdx + 1
                    && (valueIdx == -1 && rgIdx == -1 || valueIdx != rgIdx)
                    && itemCount >= (2 + (int)(scaleIdx != -1) + (int)(rgIdx != -1));
    
    if (!input_ok) {
        LOG_ERROR("ValueError", "Invalid LHA Block prototype.");
    }

    DBManager().add_lha_prototype(blockName, itemCount, valueIdx, scaleIdx, rgIdx, globalScale);
}