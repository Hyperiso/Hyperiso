#include "ParamBlockLoader.h"

static double value_to_double(const DBNode::Value& v, double fallback = 0.0)
{
    if (std::holds_alternative<double>(v))
        return std::get<double>(v);
    if (std::holds_alternative<int>(v))
        return static_cast<double>(std::get<int>(v));
    return fallback;
}

void ParamBlockLoader::load(std::shared_ptr<BlockAccessor> dest, fs::path src_file, bool block_in_blocks) {
    LOG_DEBUG("Loading parameter blocks from", src_file.string());
    auto np  = DBNodeProviderFactory::createDBNodeProvider(src_file);
    auto src = np->provide_db_as_node();
    
    if (block_in_blocks) {
        std::cout << src_file.string() << std::endl;
        for (auto &bk : src->get_keys()) {
            auto exp = src->getGroup({bk});
            for (auto &group_pair : exp) {
                auto group = src->getGroup({bk, group_pair.first});

                auto block = std::make_shared<Block>();
                block->blockname = bk + "_" + group_pair.first;
                for (auto &vk : group) {
                    const auto& key = vk.first;
                    const auto& val = vk.second;

                    if (key == "scale") {
                        if (!block->has_scale()) {
                            if (std::holds_alternative<double>(val)) {
                                block->set_scale(std::get<double>(val));
                            } else if (std::holds_alternative<int>(val)) {
                                block->set_scale(static_cast<double>(std::get<int>(val)));
                            } else {
                                LOG_WARN("ParamBlockLoader", "Non-numeric block scale under ", bk);
                            }
                        }
                        continue; 
                    }

                    if (!std::holds_alternative<std::shared_ptr<DBNode>>(val)) {
                        LOG_WARN("ParamBlockLoader", "Unexpected non-node entry under ", bk, " key ", key, " — skipping");
                        continue;
                    }

                    auto node = std::get<std::shared_ptr<DBNode>>(val);

                    if (!node->contains("central_value")) {
                        LOG_ERROR("ParamBlockLoader", "DBNode doesn't have all necessary keys (central_value).");
                        continue;
                    }

                    auto value = node->get("central_value");
                    if (std::holds_alternative<BlockName>(value)) {
                        continue;
                    }

                    double val_central = std::get<double>(value);

                    auto stat = node->contains("stat_error") ? node->get("stat_error") : DBNode::Value{0.0};
                    auto syst = node->contains("syst_error") ? node->get("syst_error") : DBNode::Value{0.0};

                    double stat_d = std::holds_alternative<double>(stat) ? std::get<double>(stat)
                                    : std::holds_alternative<int>(stat) ? static_cast<double>(std::get<int>(stat))
                                                                        : 0.0;
                    double syst_d = std::holds_alternative<double>(syst) ? std::get<double>(syst)
                                    : std::holds_alternative<int>(syst) ? static_cast<double>(std::get<int>(syst))
                                                                        : 0.0;

                    block->store(
                        LhaID(vk.first),
                        std::make_shared<Parameter>(Parameter(ParamId(bk, LhaID(vk.first)),
                                                            val_central, stat_d, syst_d))
                    );

                    if (node->contains("scale")) {
                        auto scale = node->get("scale");
                        if (std::holds_alternative<double>(scale))
                            block->retrieve(LhaID(vk.first))->set_scale(std::get<double>(scale));
                        else if (std::holds_alternative<int>(scale))
                            block->retrieve(LhaID(vk.first))->set_scale(static_cast<double>(std::get<int>(scale)));
                    }

                    if (node->contains("bin_low") || node->contains("bin_high")) {
                        if (!(node->contains("bin_low") && node->contains("bin_high")))
                            LOG_ERROR("LogicError", "Missing one end of the binning.");

                        auto bin_low = node->get("bin_low");
                        auto bin_high = node->get("bin_high");

                        double d_bin_low  = value_to_double(bin_low);
                        double d_bin_high = value_to_double(bin_high);

                        block->retrieve(LhaID(vk.first))->set_bin(std::pair(d_bin_low, d_bin_high));
                    }
                }
                dest->emplace(block->blockname, block);
            }
        }
        return;
    }

    for (auto &bk : src->get_keys()) {
        auto block = std::make_shared<Block>();
        block->blockname = bk;
        LOG_DEBUG("Loading block", bk);

        auto group = src->getGroup({bk});  

        for (auto &vk : group) {
            const auto& key = vk.first;
            const auto& val = vk.second;

            if (key == "scale") {
                if (!block->has_scale()) {
                    if (std::holds_alternative<double>(val)) {
                        block->set_scale(std::get<double>(val));
                    } else if (std::holds_alternative<int>(val)) {
                        block->set_scale(static_cast<double>(std::get<int>(val)));
                    } else {
                        LOG_WARN("ParamBlockLoader", "Non-numeric block scale under ", bk);
                    }
                }
                continue; 
            }

            if (!std::holds_alternative<std::shared_ptr<DBNode>>(val)) {
                LOG_WARN("ParamBlockLoader", "Unexpected non-node entry under ", bk, " key ", key, " — skipping");
                continue;
            }

            auto node = std::get<std::shared_ptr<DBNode>>(val);

            if (!node->contains("central_value")) {
                LOG_ERROR("ParamBlockLoader", "DBNode doesn't have all necessary keys (central_value).");
                continue;
            }

            auto value = node->get("central_value");
            if (std::holds_alternative<BlockName>(value)) {
                continue;
            }

            double val_central = std::get<double>(value);

            auto stat = node->contains("stat_error") ? node->get("stat_error") : DBNode::Value{0.0};
            auto syst = node->contains("syst_error") ? node->get("syst_error") : DBNode::Value{0.0};

            double stat_d = std::holds_alternative<double>(stat) ? std::get<double>(stat)
                             : std::holds_alternative<int>(stat) ? static_cast<double>(std::get<int>(stat))
                                                                 : 0.0;
            double syst_d = std::holds_alternative<double>(syst) ? std::get<double>(syst)
                             : std::holds_alternative<int>(syst) ? static_cast<double>(std::get<int>(syst))
                                                                 : 0.0;

            block->store(
                LhaID(vk.first),
                std::make_shared<Parameter>(Parameter(ParamId(bk, LhaID(vk.first)),
                                                      val_central, stat_d, syst_d))
            );

            if (node->contains("scale")) {
                auto scale = node->get("scale");
                if (std::holds_alternative<double>(scale))
                    block->retrieve(LhaID(vk.first))->set_scale(std::get<double>(scale));
                else if (std::holds_alternative<int>(scale))
                    block->retrieve(LhaID(vk.first))->set_scale(static_cast<double>(std::get<int>(scale)));
            }

            if (node->contains("bin_low") || node->contains("bin_high")) {
                if (!(node->contains("bin_low") && node->contains("bin_high")))
                    LOG_ERROR("LogicError", "Missing one end of the binning.");

                auto bin_low = node->get("bin_low");
                auto bin_high = node->get("bin_high");

                double d_bin_low  = value_to_double(bin_low);
                double d_bin_high = value_to_double(bin_high);

                block->retrieve(LhaID(vk.first))->set_bin(std::pair(d_bin_low, d_bin_high));
            }
        }

        dest->emplace(bk, block);
    }

    LOG_DEBUG("Parameter blocks loaded");
}
