#pragma once
#include "WilsonGroup.h"

class GenericCoefficientGroup : public CoefficientGroup {
public:
    using CoefficientGroup::CoefficientGroup;

    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<GenericCoefficientGroup>(*this);
    }

    std::shared_ptr<CoefficientGroup> get_sm_group() override {
        auto g = std::make_shared<GenericCoefficientGroup>(*this);
        g->set_wilson_type(ContributionType::SM);
        return g;
    }
};
