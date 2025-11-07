#ifndef CUSTOM_WILSON_H
#define CUSTOM_WILSON_H

#include "Wilson.h"

class CustomWilson : public WilsonCoefficient {
public:

    CustomWilson(const LhaID& name,
                 const std::string& storage_block,
                 QCDOrder max_order = QCDOrder::LO,
                 ContributionType type = ContributionType::SM)
    : WilsonCoefficient(name, storage_block, type)
    {
        std::cout << "Creating CustomWilson: " << name << std::endl;
        matching_info[QCDOrder::LO] = MatchingInfo();
        matching_info[QCDOrder::NLO] = MatchingInfo();
        matching_info[QCDOrder::NNLO] = MatchingInfo();

        // this->max_order = max_order;
        this->type = type;
    }

    CustomWilson& set_order_info(
        QCDOrder order,
        std::unordered_set<ParamId> sources,
        std::function<scalar_t(const ParamSrc&)> compute,
        LhaID lhaid
    ) {
        matching_info[order] = MatchingInfo(std::move(sources), std::move(compute), std::move(lhaid));
        // if (order > max_order) max_order = order;
        return *this;
    }

    CustomWilson& with_type(ContributionType t) { this->type = t; return *this; }

    CustomWilson& with_storage_block(std::string block) { this->storage_block = std::move(block); return *this; }

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CustomWilson>(*this);
    }
};

#endif // CUSTOM_WILSON_H
