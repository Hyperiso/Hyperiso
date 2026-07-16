#ifndef USER_PARAMETER_PROXY_H
#define USER_PARAMETER_PROXY_H

#include "IUserParameterProxy.h"
#include "ParameterProvider.h"
#include "ParameterSetter.h"
#include "Include.h"

class UserParameterProxy : public IUserParameterProxy<std::string, LhaID> {
public:

    UserParameterProxy(std::vector<ParameterType> types);


    std::optional<double> get_value(std::string block, LhaID id) override;

    void set_value(std::string block, LhaID id, double val) override;

private:
    std::optional<ParameterType> get_proxy(std::string block, LhaID id);

    std::map<ParameterType, ParameterProvider> types_p;

};

#endif 