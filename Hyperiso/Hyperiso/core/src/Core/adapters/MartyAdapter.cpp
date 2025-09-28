#include "MartyAdapter.h"

fs::path MartyAdapter::get_path(MartyPath path_name) {

    std::shared_ptr<DefaultPathsProvider> dpp = std::make_shared<DefaultPathsProvider>();

    switch (path_name) {
    case MartyPath::MODEL_FILE:
        return MemoryManager::GetInstance()->getMemoryCache().config.mty_model_path.value();
        break;
    case MartyPath::PARAM_MAPPING_DIR:
        return dpp->param_mapping_dir_path();
        break;
    case MartyPath::TEMPLATE_DIR:
        return dpp->template_dir_path();
        break;
    default:
        LOG_ERROR("ValueError", "Unknown path for MartyAdapter.");
    };
}

bool MartyAdapter::check_flag(InternalFlag flag) {
    return MemoryManager::GetInstance()->getMemoryCache().flags.at(flag);
}

std::string MartyAdapter::get_marty_model_name() const{
    auto m_name = MemoryManager::GetInstance()->getMemoryCache().config.mty_model_name;

    if (!m_name.has_value()) {
        LOG_ERROR("LogicError", "MARTY Model name is not defined.");
    }

    return m_name.value();
}
