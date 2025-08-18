#ifndef CUSTOM_WILSON_GROUP_H
#define CUSTOM_WILSON_GROUP_H

#include "WilsonGroup.h"

class CustomCoefficientGroup : public CoefficientGroup {
public:

    explicit CustomCoefficientGroup(WGroup id,
                                    ContributionType type = ContributionType::SM)
    {
        this->id = id;
        this->wilson_type = type;
    }

    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<CustomCoefficientGroup>(*this);
    }

    CustomCoefficientGroup& add_coefficient(const std::shared_ptr<WilsonCoefficient>& coef) {
        this->insert({coef->get_name(), coef});
        return *this;
    }

    CustomCoefficientGroup& set_basis_order_sources_and_running(
        WilsonBasis basis,
        QCDOrder order,
        std::unordered_map<ParameterType, std::vector<std::string>> source_names,
        std::function<std::unordered_map<WCoef, scalar_t>(
            const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>&,
            const std::unordered_map<std::string, std::shared_ptr<Block>>&
        )> running_func
    ) {
        sources[basis][order].sources = std::move(source_names);
        sources[basis][order].func    = std::move(running_func);
        if (order > current_order) current_order = order;
        return *this;
    }

    static std::unordered_map<WCoef, scalar_t> identity_running(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
        const std::unordered_map<std::string, std::shared_ptr<Block>>& /*src*/
    ) {
        QCDOrder best = QCDOrder::LO;
        if (coef_matching.count(QCDOrder::NNLO)) best = QCDOrder::NNLO;
        else if (coef_matching.count(QCDOrder::NLO)) best = QCDOrder::NLO;

        return coef_matching.at(best);
    }

    CustomCoefficientGroup& finalize(QCDOrder order) {
        this->init(order);
        return *this;
    }
};

#endif // CUSTOM_WILSON_GROUP_H
