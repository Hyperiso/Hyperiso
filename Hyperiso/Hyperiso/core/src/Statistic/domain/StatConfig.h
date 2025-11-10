#include "AbstractConfig.h"
#include "General.h"

class StatConfig : public AbstractConfig {
    std::vector<ParamId> nuisances;
    std::vector<ParamId> parameters;
};