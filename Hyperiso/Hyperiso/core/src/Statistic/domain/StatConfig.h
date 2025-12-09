#include "AbstractConfig.h"
#include "ParamID.h"

class StatConfig : public AbstractConfig {
    std::vector<ParamId> nuisances;
    std::vector<ParamId> parameters;
};