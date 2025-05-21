#ifndef WILSONINTERFACE_H
#define WILSONINTERFACE_H

#include "WilsonManager.h"
#include "Wilson_SUSY.h"
#include "Wilson_THDM_super.h"
#include "Wilson_SUSY_super.h"
#include "MartyWilsonSuper.h"
#include "General.h"
#include "AbstractConfig.h"

#include <map>

class WilsonInterface {
private:
    CoefficientManager wm;
    std::map<std::string, std::shared_ptr<CoefficientGroup>> group_ptrs;

    QCDOrder ensure_mty_compat(QCDOrder order) {
        if (UseMarty().get() && !(order == QCDOrder::LO)) {
            LOG_WARN("Using MARTY defaults all calculations to LO in QCD.");
            return QCDOrder::LO;
        }
        return order;
    }

public:
    explicit WilsonInterface() {
        LOG_INFO("In WilsonInterface constructor");
        this->init();
        this->wm = CoefficientManager();
    }

    void init() {
        LOG_INFO("In WilsonInterface::init");
        WilsonParameterHelper().init(2);

        this->group_ptrs = { };

        std::map<std::string, std::shared_ptr<CoefficientGroup>> bsm_groups;
        if (ModelAPI().get() == Model::SM) {
            bsm_groups = {
                {"BCoefficients", std::make_shared<BCoefficientGroup>()},
                {"BPrimeCoefficients", std::make_shared<BPrimeCoefficientGroup>()},
                {"BScalarCoefficients", std::make_shared<BScalarCoefficientGroup>()},
                {"BlnuCoefficients", std::make_shared<BlnuCoefficientGroup>()},
                {"BclnuCoefficients", std::make_shared<BclnuCoefficientGroup>()}
            };
        } else if (ModelAPI().get() == Model::THDM) {
            bsm_groups = {
                {"BCoefficients", std::make_shared<BCoefficientGroup_THDM>()},
                {"BPrimeCoefficients", std::make_shared<BPrimeCoefficientGroup_THDM>()},
                {"BScalarCoefficients", std::make_shared<BScalarCoefficientGroup_THDM>()},
                {"BlnuCoefficients", std::make_shared<BlnuCoefficientGroup_THDM>()},
                {"BclnuCoefficients", std::make_shared<BclnuCoefficientGroup_THDM>()}
            };
        } else if (ModelAPI().get() == Model::SUSY) {
            bsm_groups = {
                {"BCoefficients", std::make_shared<BCoefficientGroup_susy>()},
                {"BPrimeCoefficients", std::make_shared<BPrimeCoefficientGroup_susy>()},
                {"BScalarCoefficients", std::make_shared<BScalarCoefficientGroup_susy>()},
                {"BlnuCoefficients", std::make_shared<BlnuCoefficientGroup_SUSY>()},
                {"BclnuCoefficients", std::make_shared<BclnuCoefficientGroup_SUSY>()}
            };
        }

        this->group_ptrs.insert(bsm_groups.begin(), bsm_groups.end());
    }

    void addWilsonGroup(WGroup group_name) {
        this->wm.registerCoefficientGroup(GroupMapper::str(group_name), this->group_ptrs.at(GroupMapper::str(group_name)));
    }

    void set_matching_scale(double mu_W) {
        ScaleSetter(ScaleType::MATCHING).set(mu_W);
    }

    void init_group_matching(WGroup group_name, QCDOrder order) {
        this->wm.init_group_matching(GroupMapper::str(group_name), OrderMapper::str(ensure_mty_compat(order)));
    }

    void init_group_hadronic(WGroup group_name, QCDOrder order) {
        this->wm.init_group_hadronic(GroupMapper::str(group_name), OrderMapper::str(ensure_mty_compat(order)));
    }

    void set_hadronic_scale(double mu_h) {
        ScaleSetter(ScaleType::HADRONIC).set(mu_h);
    }

    void switchbasis(WGroup group_name) {
        this->wm.switchbasis(GroupMapper::str(group_name));
    }

    complex_t getMatchingCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
        return this->wm.getMatchingCoefficient(GroupMapper::str(group), WCoefMapper::str(coeff), OrderMapper::str(order), cont_type);
    }

    complex_t getM(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
        return getMatchingCoefficient(group, coeff, order, cont_type);
    }

    complex_t getFullMatchingCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
        return this->wm.getFullMatchingCoefficient(GroupMapper::str(group), WCoefMapper::str(coeff), OrderMapper::str(ensure_mty_compat(order)), cont_type);
    }

    complex_t getFM(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
        return getFullMatchingCoefficient(group, coeff, order, cont_type);
    }

    complex_t getRunCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
        return this->wm.getRunCoefficient(GroupMapper::str(group), WCoefMapper::str(coeff), OrderMapper::str(order), cont_type);
    }

    complex_t getR(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
        return getRunCoefficient(group, coeff, order, cont_type);
    }

    complex_t getFullRunCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
        return this->wm.getFullRunCoefficient(GroupMapper::str(group), WCoefMapper::str(coeff), OrderMapper::str(ensure_mty_compat(order)), cont_type);
    }

    complex_t getFR(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
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

    std::map<QCDOrder, complex_t> getSepOrderRunCoefficient(WGroup group, WCoef coeff, ContributionType cont_type) {
        std::map<QCDOrder, complex_t> C {{
            {QCDOrder::LO, getRunCoefficient(group, coeff, QCDOrder::LO, cont_type)},
            {QCDOrder::NLO, getRunCoefficient(group, coeff, QCDOrder::NLO, cont_type)},
            {QCDOrder::NNLO, getRunCoefficient(group, coeff, QCDOrder::NNLO, cont_type)}
        }};
        return C;
    }

    std::map<QCDOrder, complex_t> getSR(WGroup group, WCoef coeff, ContributionType cont_type) {
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

    std::map<WCoef, complex_t> getAllRunCoefficients(WGroup group, QCDOrder order, ContributionType cont_type) {
        std::vector<WCoef> ids = WCoefMapper::get_group(group);
        std::map<WCoef, complex_t> Cs;
        for (auto c : ids) {
            Cs.emplace(c, getRunCoefficient(group, c, order, cont_type));
        }

        return Cs;
    }

    std::map<WCoef, complex_t> getAR(WGroup group, QCDOrder order, ContributionType cont_type) {
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

    std::map<WCoef, complex_t> getAllFullRunCoefficients(WGroup group, QCDOrder order, ContributionType cont_type) {
        std::vector<WCoef> ids = WCoefMapper::get_group(group);
        std::map<WCoef, complex_t> Cs;
        for (auto c : ids) {
            Cs.emplace(c, getFullRunCoefficient(group, c, order, cont_type));
        }

        return Cs;
    }

    std::map<WCoef, complex_t> getAFR(WGroup group, QCDOrder order, ContributionType cont_type) {
        return getAllFullRunCoefficients(group, order, cont_type);
    }
};

#endif // __WILSONINTERFACE_H__