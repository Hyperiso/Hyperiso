#include "LhaNodeProvider.h"

static bool valid_opt_idx(int idx, size_t itemCount) {
    return idx == -1 || (idx >= 0 && static_cast<size_t>(idx) < itemCount);
}

static bool conflicts_with_value(int idx, size_t valueIdx) {
    return idx >= 0 && static_cast<size_t>(idx) == valueIdx;
}

static bool conflicts(int a, int b) {
    return a >= 0 && b >= 0 && a == b;
}

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
    bool input_ok = true;

    input_ok &= (valueIdx < itemCount);

    input_ok &= valid_opt_idx(scaleIdx, itemCount);
    input_ok &= valid_opt_idx(rgIdx, itemCount);
    input_ok &= valid_opt_idx(binIdx, itemCount);

    if (binIdx != -1) {
        input_ok &= (binIdx >= 0);
        input_ok &= (static_cast<size_t>(binIdx + 1) < itemCount);
    }

    input_ok &= !conflicts_with_value(scaleIdx, valueIdx);
    input_ok &= !conflicts_with_value(rgIdx, valueIdx);
    input_ok &= !conflicts_with_value(binIdx, valueIdx);
    if (binIdx != -1) {
        input_ok &= (static_cast<size_t>(binIdx + 1) != valueIdx);
    }

    input_ok &= !conflicts(scaleIdx, rgIdx);
    input_ok &= !conflicts(scaleIdx, binIdx);
    input_ok &= !(binIdx != -1 && scaleIdx >= 0 && scaleIdx == binIdx + 1);
    input_ok &= !conflicts(rgIdx, binIdx);
    input_ok &= !(binIdx != -1 && rgIdx >= 0 && rgIdx == binIdx + 1);

    const size_t minCols =
        2
        + (scaleIdx != -1 ? 1 : 0)
        + (rgIdx    != -1 ? 1 : 0)
        + (binIdx   != -1 ? 2 : 0);

    input_ok &= (itemCount >= minCols);

    if (!input_ok) {
        LOG_ERROR("ValueError", "Invalid LHA Block prototype.");
        return;
    }

    DBManager().add_lha_prototype(blockName, itemCount, valueIdx, scaleIdx, rgIdx, binIdx, globalScale);
}