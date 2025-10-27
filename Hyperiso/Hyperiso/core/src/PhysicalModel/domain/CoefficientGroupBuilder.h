#pragma once
#include "GroupDefinition.h"
#include "WilsonGroup.h"
#include "Wilson.h"
#include "WilsonCoefficientRegistry.h"
#include "GenericWilsonGroup.h"

class CoefficientGroupBuilder {
public:
    explicit CoefficientGroupBuilder(const CoefficientRegistry& reg) : reg_(reg) {}

    std::shared_ptr<CoefficientGroup> build(const BuildContext& ctx) const;

private:
    const CoefficientRegistry& reg_;
};
