#include "ConfigProvider.h"

bool ConfigProvider::check_flag(ExternalFlag flag) {
    return MemoryManager::GetInstance()->getMemoryCache().config.flags.at(flag);
}

Model ConfigProvider::get_model() {
    return MemoryManager::GetInstance()->getMemoryCache().config.model;
}
