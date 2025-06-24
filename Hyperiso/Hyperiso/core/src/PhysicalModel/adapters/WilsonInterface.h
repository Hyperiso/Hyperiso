#ifndef WILSONINTERFACE_H
#define WILSONINTERFACE_H

#include "WilsonManager.h"
#include "Wilson_SUSY.h"
#include "Wilson_THDM_super.h"
#include "Wilson_SUSY_super.h"
#include "MartyWilsonSuper.h"
#include "General.h"
#include "AbstractConfig.h"
#include "WilsonBuilder.h"

#include <map>

class WilsonInterface {
private:
    std::shared_ptr<WilsonBuilder> builder;
    std::shared_ptr<WilsonProvider> provider;

    QCDOrder ensure_mty_compat(QCDOrder order) {
        if (UseMarty().get() && !(order == QCDOrder::LO)) {
            LOG_WARN("Using MARTY defaults all calculations to LO in QCD.");
            return QCDOrder::LO;
        }
        return order;
    }

    bool built {false};

public:
    WilsonInterface() = default;

    void build(WilsonBuildConfig config) {
        this->builder = std::make_shared<WilsonBuilder>(config);
        this->provider = this->builder->get_wilson_provider();
        built = true;
    }

    // TODO
    // void addWilsonGroup(WGroup group_id) {
    //     this->builder->add(GroupMapper::str(group_id), this->group_ptrs.at(GroupMapper::str(group_name)));
    // }

    void set_matching_scale(double mu_W) {
        ScaleSetter(ScaleType::MATCHING).set(mu_W);
    }

    void set_hadronic_scale(double mu_h) {
        ScaleSetter(ScaleType::HADRONIC).set(mu_h);
    }

    complex_t getMatchingCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
        if (!built) {
            LOG_ERROR("LogicError", "Interface has not been built");
        }

        WilsonRequest request {
            group,
            coeff,
            order,
            cont_type,
            ScaleType::MATCHING,
            false // sum_qcd_orders
        };
        return this->provider->get(std::make_shared<WilsonRequest>(request));
    }

    complex_t getM(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
        return getMatchingCoefficient(group, coeff, order, cont_type);
    }

    complex_t getFullMatchingCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
        if (!built) {
            LOG_ERROR("LogicError", "Interface has not been built");
        }

        WilsonRequest request {
            group,
            coeff,
            order,
            cont_type,
            ScaleType::MATCHING,
            true // sum_qcd_orders
        };
        return this->provider->get(std::make_shared<WilsonRequest>(request));
    }

    complex_t getFM(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
        return getFullMatchingCoefficient(group, coeff, order, cont_type);
    }

    complex_t getRunCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD) {
        if (!built) {
            LOG_ERROR("LogicError", "Interface has not been built");
        }

        WilsonRequest request {
            group,
            coeff,
            order,
            cont_type,
            ScaleType::HADRONIC,
            false // sum_qcd_orders
        };
        request.basis = basis;
        return this->provider->get(std::make_shared<WilsonRequest>(request));
    }

    complex_t getR(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD) {
        return getRunCoefficient(group, coeff, order, cont_type);
    }

    complex_t getFullRunCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD) {
        if (!built) {
            LOG_ERROR("LogicError", "Interface has not been built");
        }

        WilsonRequest request {
            group,
            coeff,
            order,
            cont_type,
            ScaleType::HADRONIC,
            true // sum_qcd_orders
        };
        request.basis = basis;
        return this->provider->get(std::make_shared<WilsonRequest>(request));
    }

    complex_t getFR(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD) {
        return getFullRunCoefficient(group, coeff, order, cont_type);
    }

    std::map<QCDOrder, complex_t> getSepOrderMatchingCoefficient(WGroup group, WCoef coeff, ContributionType cont_type) {
        std::map<QCDOrder, complex_t> C {{
            {QCDOrder::LO, getMatchingCoefficient(group, coeff, QCDOrder::LO, cont_type)},
            {QCDOrder::NLO, getMatchingCoefficient(group, coeff, QCDOrder::NLO, cont_type)},
            {QCDOrder::NNLO, getMatchingCoefficient(group, coeff, QCDOrder::NNLO, cont_type)}
        }};
        return C;
    }

    std::map<QCDOrder, complex_t> getSM(WGroup group, WCoef coeff, ContributionType cont_type) {
        return getSepOrderMatchingCoefficient(group, coeff, cont_type);
    }

    std::map<QCDOrder, complex_t> getSepOrderRunCoefficient(WGroup group, WCoef coeff, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD) {
        std::map<QCDOrder, complex_t> C {{
            {QCDOrder::LO, getRunCoefficient(group, coeff, QCDOrder::LO, cont_type)},
            {QCDOrder::NLO, getRunCoefficient(group, coeff, QCDOrder::NLO, cont_type)},
            {QCDOrder::NNLO, getRunCoefficient(group, coeff, QCDOrder::NNLO, cont_type)}
        }};
        return C;
    }

    std::map<QCDOrder, complex_t> getSR(WGroup group, WCoef coeff, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD) {
        return getSepOrderRunCoefficient(group, coeff, cont_type);
    }

    std::map<WCoef, complex_t> getAllMatchingCoefficients(WGroup group, QCDOrder order, ContributionType cont_type) {
        std::vector<WCoef> ids = WCoefMapper::get_group(group);
        std::map<WCoef, complex_t> Cs;
        for (auto c : ids) {
            Cs.emplace(c, getMatchingCoefficient(group, c, order, cont_type));
        }

        return Cs;
    }

    std::map<WCoef, complex_t> getAM(WGroup group, QCDOrder order, ContributionType cont_type) {
        return getAllMatchingCoefficients(group, order, cont_type);
    }

    std::map<WCoef, complex_t> getAllRunCoefficients(WGroup group, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD) {
        std::vector<WCoef> ids = WCoefMapper::get_group(group);
        std::map<WCoef, complex_t> Cs;
        for (auto c : ids) {
            Cs.emplace(c, getRunCoefficient(group, c, order, cont_type));
        }

        return Cs;
    }

    std::map<WCoef, complex_t> getAR(WGroup group, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD) {
        return getAllRunCoefficients(group, order, cont_type);
    }

    std::map<WCoef, complex_t> getAllFullMatchingCoefficients(WGroup group, QCDOrder order, ContributionType cont_type) {
        std::vector<WCoef> ids = WCoefMapper::get_group(group);
        std::map<WCoef, complex_t> Cs;
        for (auto c : ids) {
            Cs.emplace(c, getFullMatchingCoefficient(group, c, order, cont_type));
        }

        return Cs;
    }

    std::map<WCoef, complex_t> getAFM(WGroup group, QCDOrder order, ContributionType cont_type) {
        return getAllFullMatchingCoefficients(group, order, cont_type);
    }

    std::map<WCoef, complex_t> getAllFullRunCoefficients(WGroup group, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD) {
        std::vector<WCoef> ids = WCoefMapper::get_group(group);
        std::map<WCoef, complex_t> Cs;
        for (auto c : ids) {
            Cs.emplace(c, getFullRunCoefficient(group, c, order, cont_type));
        }

        return Cs;
    }

    std::map<WCoef, complex_t> getAFR(WGroup group, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD) {
        return getAllFullRunCoefficients(group, order, cont_type);
    }
};

#endif // __WILSONINTERFACE_H__