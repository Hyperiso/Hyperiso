#include "MartyAdapter.h"

#include <stdexcept>

fs::path MartyAdapter::get_path(MartyPath path_name) {
    auto* mm = MemoryManager::GetInstance();

    switch (path_name) {
    case MartyPath::MODEL_FILE: {
        const auto& config = mm->getMemoryCache().config;
        if (!config.mty_model_path.has_value()) {
            throw std::logic_error("MARTY model path is not defined in HyperisoConfig::mty_model_path.");
        }
        return config.mty_model_path.value();
    }

    case MartyPath::PARAM_MAPPING_DIR:
        return mm->get_path(APIPath::PARAM_MAPPING_DIR);

    case MartyPath::SM_MAPPING_FILE:
        return mm->get_path(APIPath::PARAM_MAPPING_DIR) / "sm.json";

    case MartyPath::TEMPLATE_DIR:
        return mm->get_path(APIPath::TEMPLATE_DIR) / "MARTY";

    case MartyPath::MARTY_TEMP_DIR:
        return mm->get_path(APIPath::MARTY_TEMP_DIR);

    case MartyPath::BSM_MAPPING_FILE: {
        auto optional_path = get_optional_path(path_name);
        if (!optional_path.has_value()) {
            throw std::logic_error("MARTY BSM mapping path is not defined in HyperisoConfig::mty_bsm_mapping_path.");
        }
        return optional_path.value();
    }
    }

    throw std::logic_error("Unknown MartyPath.");
}

std::optional<fs::path> MartyAdapter::get_optional_path(MartyPath path_name) {
    if (path_name != MartyPath::BSM_MAPPING_FILE) {
        return get_path(path_name);
    }

    const auto& config = MemoryManager::GetInstance()->getMemoryCache().config;
    return config.mty_bsm_mapping_path;
}

bool MartyAdapter::check_flag(InternalFlag flag) {
    return MemoryManager::GetInstance()->getMemoryCache().flags.at(flag);
}

std::string MartyAdapter::get_marty_model_name() const {
    const auto& config = MemoryManager::GetInstance()->getMemoryCache().config;

    if (!config.mty_model_name.has_value()) {
        throw std::logic_error("MARTY model name is not defined in HyperisoConfig::mty_model_name.");
    }

    return config.mty_model_name.value();
}
