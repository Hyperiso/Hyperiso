#include "UserParameterProxy.h"


UserParameterProxy::UserParameterProxy(std::vector<ParameterType> types) {
    for (auto t : types) {
        this->types_p[t] = ParameterProvider(t);
    }
}
std::optional<double> UserParameterProxy::get_value(std::string block, LhaID id) {
    std::optional<ParameterType> pt = this->get_proxy(block, id);

    if (!pt.has_value()) {
        return std::nullopt;
    }

    return std::optional<double>(this->types_p[pt.value()](block, id));
}

void UserParameterProxy::set_value(std::string block, LhaID id, double val) {
    std::optional<ParameterType> pt = this->get_proxy(block, id);

    if (!pt.has_value()) {
        LOG_ERROR("ValueError", "Cannot set a parameter :", block, id, "that does not exist.");
    }

    ParameterSetter().mutate({pt.value(), block, id}, val);
}

std::optional<ParameterType> UserParameterProxy::get_proxy(std::string block, LhaID id) {
    for (auto p : this->types_p) {
        if (p.second.exists(block, id)) {
            return std::optional<ParameterType>(p.first);
        }
    }
    return std::nullopt;
}