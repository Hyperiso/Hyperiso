#include "ParamBlockWriter.h"

static BlockName lhaid_to_key(const LhaID& id)
{
    std::ostringstream oss;
    oss << id;
    return BlockName(oss.str());
}

void ParamBlockWriter::write(std::shared_ptr<DBNode> dest, std::shared_ptr<BlockAccessor> src) {

    if (!dest) {
        LOG_ERROR("ParamBlockWriter", "dest is null");
        return;
    }
    if (!src) {
        LOG_ERROR("ParamBlockWriter", "src is null");
        return;
    }

    for (const auto& bk : src->get_block_names()) {

        if (!src->contains(bk)) {
            LOG_WARN("ParamBlockWriter", "Block not found in accessor:", bk);
            continue;
        }

        auto blk = src->at(bk);
        if (!blk) {
            LOG_WARN("ParamBlockWriter", "Null block pointer for:", bk);
            continue;
        }

        std::map<BlockName, DBNode::Value> groupData;

        if (blk->has_scale()) {
            groupData[BlockName("scale")] = blk->get_scale();
        }

        for (const auto& [id, p] : blk->getItems()) {
            if (!p) continue;

            auto node = std::make_shared<DBNode>();

            const scalar_t value = p->get_val();
            node->set(value.real(), "central_value");
            if (value.imag() != 0.0) {
                node->set(value.imag(), "imaginary_value");
            }

            {
                auto [stat, syst] = p->get_std();
                node->set(static_cast<double>(stat), "stat_error");
                node->set(static_cast<double>(syst), "syst_error");
            }

            {
                double sc = p->get_scale();
                if (sc != -1.0) {
                    node->set(sc, "scale");
                }
            }

            {
                auto bin = p->get_bin(); 
                if (!(bin.first == -1.0 && bin.second == -1.0)) {
                    node->set(bin.first,  "bin_low");
                    node->set(bin.second, "bin_high");
                }
            }

            groupData[lhaid_to_key(id)] = node;
        }

        dest->setGroup({bk}, groupData);
    }
}
