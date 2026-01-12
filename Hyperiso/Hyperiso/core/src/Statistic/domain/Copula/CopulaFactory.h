#ifndef __COPULAFACTORY_H__
#define __COPULAFACTORY_H__

#include <memory>
#include <random>
#include <algorithm>
#include <variant>

#include "CopulaType.h"
#include "ICopula.h"
#include "GaussianCopula.h"
#include "StudentTCopula.h"

using CopulaConfig = std::variant<GaussianCopulaConfig, StudentTCopulaConfig>;

class CopulaFactory {
public:
    static std::unique_ptr<ICopula> create(CopulaType name, CopulaConfig config, unsigned int seed = std::random_device{}());
};

#endif // __COPULAFACTORY_H__
