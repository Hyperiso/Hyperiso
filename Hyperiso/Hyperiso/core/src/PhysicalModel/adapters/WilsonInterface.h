#ifndef WILSONINTERFACE_H
#define WILSONINTERFACE_H

#include "WilsonManager.h"
#include "Wilson_SUSY.h"
#include "Wilson_THDM.h"
#include "MartyWilson.h"
#include "General.h"

#include <map>

struct WilsonConfig {
    std::unordered_set<WGroup> groups;
    double matching_scale;
    double hadronic_scale;
    QCDOrder order{QCDOrder::LO}; //TODO : I choose LO as default here
};

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
        this->init();
        this->wm = CoefficientManager();
    }

    void init() {
        WilsonParameterHelper().init(2);

        this->group_ptrs = {
            {"BCoefficients", std::make_shared<BCoefficientGroup>()},
            {"BPrimeCoefficients", std::make_shared<BPrimeCoefficientGroup>()},
            {"BScalarCoefficients", std::make_shared<BScalarCoefficientGroup>()},
            {"BlnuCoefficients", std::make_shared<BlnuCoefficientGroup>()},
            {"BclnuCoefficients", std::make_shared<BclnuCoefficientGroup>()}
        };

        std::map<std::string, std::shared_ptr<CoefficientGroup>> bsm_groups;
        if (ModelAPI().get() == Model::THDM) {
            bsm_groups = {
                {"BCoefficients_THDM", std::make_shared<BCoefficientGroup_THDM>()},
                {"BPrimeCoefficients_THDM", std::make_shared<BPrimeCoefficientGroup_THDM>()},
                {"BScalarCoefficients_THDM", std::make_shared<BScalarCoefficientGroup_THDM>()},
                {"BlnuCoefficients_THDM", std::make_shared<BlnuCoefficientGroup_THDM>()},
                {"BclnuCoefficients_THDM", std::make_shared<BclnuCoefficientGroup_THDM>()}
            };
        } else if (ModelAPI().get() == Model::SUSY) {
            bsm_groups = {
                {"BCoefficients_SUSY", std::make_shared<BCoefficientGroup_susy>()},
                {"BPrimeCoefficients_SUSY", std::make_shared<BPrimeCoefficientGroup_susy>()},
                {"BScalarCoefficients_SUSY", std::make_shared<BScalarCoefficientGroup_susy>()},
                {"BlnuCoefficients_SUSY", std::make_shared<BlnuCoefficientGroup_SUSY>()},
                {"BclnuCoefficients_SUSY", std::make_shared<BclnuCoefficientGroup_SUSY>()}
            };
        }

        this->group_ptrs.insert(bsm_groups.begin(), bsm_groups.end());
    }

    void addWilsonGroup(WGroup group_name) {
        this->wm.registerCoefficientGroup(GroupMapper::str(group_name), this->group_ptrs.at(GroupMapper::str(group_name)));
    }

    void set_matching_scale(double mu_W) {
        this->wm.set_matching_scale(mu_W);
    }

    void init_group_matching(WGroup group_name, QCDOrder order) {
        this->wm.init_group_matching(GroupMapper::str(group_name), OrderMapper::str(ensure_mty_compat(order)));
    }

    void init_group_hadronic(WGroup group_name, QCDOrder order) {
        this->wm.init_group_hadronic(GroupMapper::str(group_name), OrderMapper::str(ensure_mty_compat(order)));
    }

    void set_hadronic_scale(WGroup group_name, double mu_h) {
        this->wm.set_hadronic_scale(mu_h);
    }

    void switchbasis(WGroup group_name) {
        this->wm.switchbasis(GroupMapper::str(group_name));
    }

    complex_t getMatchingCoefficient(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return this->wm.getMatchingCoefficient(GroupMapper::str(group), WCoefMapper::str(coeff), OrderMapper::str(order), sm_only);
    }

    complex_t getM(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return getMatchingCoefficient(group, coeff, order, sm_only);
    }

    complex_t getFullMatchingCoefficient(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return this->wm.getFullMatchingCoefficient(GroupMapper::str(group), WCoefMapper::str(coeff), OrderMapper::str(ensure_mty_compat(order)), sm_only);
    }

    complex_t getFM(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return getFullMatchingCoefficient(group, coeff, order, sm_only);
    }

    complex_t getRunCoefficient(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return this->wm.getRunCoefficient(GroupMapper::str(group), WCoefMapper::str(coeff), OrderMapper::str(order), sm_only);
    }

    complex_t getR(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return getRunCoefficient(group, coeff, order, sm_only);
    }

    complex_t getFullRunCoefficient(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return this->wm.getFullRunCoefficient(GroupMapper::str(group), WCoefMapper::str(coeff), OrderMapper::str(ensure_mty_compat(order)), sm_only);
    }

    complex_t getFR(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return getFullRunCoefficient(group, coeff, order, sm_only);
    }

    std::map<QCDOrder, complex_t> getSepOrderMatchingCoefficient(WGroup group, WCoef coeff, bool sm_only=false) {
        std::map<QCDOrder, complex_t> C {{
            {QCDOrder::LO, getMatchingCoefficient(group, coeff, QCDOrder::LO)},
            {QCDOrder::NLO, getMatchingCoefficient(group, coeff, QCDOrder::NLO)},
            {QCDOrder::NNLO, getMatchingCoefficient(group, coeff, QCDOrder::NNLO)}
        }};
        return C;
    }

    std::map<QCDOrder, complex_t> getSM(WGroup group, WCoef coeff, bool sm_only=false) {
        return getSepOrderMatchingCoefficient(group, coeff, sm_only);
    }

    std::map<QCDOrder, complex_t> getSepOrderRunCoefficient(WGroup group, WCoef coeff, bool sm_only=false) {
        std::map<QCDOrder, complex_t> C {{
            {QCDOrder::LO, getRunCoefficient(group, coeff, QCDOrder::LO, sm_only)},
            {QCDOrder::NLO, getRunCoefficient(group, coeff, QCDOrder::NLO, sm_only)},
            {QCDOrder::NNLO, getRunCoefficient(group, coeff, QCDOrder::NNLO, sm_only)}
        }};
        return C;
    }

    std::map<QCDOrder, complex_t> getSR(WGroup group, WCoef coeff, bool sm_only=false) {
        return getSepOrderRunCoefficient(group, coeff, sm_only);
    }

    std::map<WCoef, complex_t> getAllMatchingCoefficients(WGroup group, QCDOrder order, bool sm_only=false) {
        std::vector<WCoef> ids = WCoefMapper::get_group(group);
        std::map<WCoef, complex_t> Cs;
        for (auto c : ids) {
            Cs.emplace(c, getMatchingCoefficient(group, c, order, sm_only));
        }

        return Cs;
    }

    std::map<WCoef, complex_t> getAM(WGroup group, QCDOrder order, bool sm_only=false) {
        return getAllMatchingCoefficients(group, order, sm_only);
    }

    std::map<WCoef, complex_t> getAllRunCoefficients(WGroup group, QCDOrder order, bool sm_only=false) {
        std::vector<WCoef> ids = WCoefMapper::get_group(group);
        std::map<WCoef, complex_t> Cs;
        for (auto c : ids) {
            Cs.emplace(c, getRunCoefficient(group, c, order, sm_only));
        }

        return Cs;
    }

    std::map<WCoef, complex_t> getAR(WGroup group, QCDOrder order, bool sm_only=false) {
        return getAllRunCoefficients(group, order, sm_only);
    }

    std::map<WCoef, complex_t> getAllFullMatchingCoefficients(WGroup group, QCDOrder order, bool sm_only=false) {
        std::vector<WCoef> ids = WCoefMapper::get_group(group);
        std::map<WCoef, complex_t> Cs;
        for (auto c : ids) {
            Cs.emplace(c, getFullMatchingCoefficient(group, c, order, sm_only));
        }

        return Cs;
    }

    std::map<WCoef, complex_t> getAFM(WGroup group, QCDOrder order, bool sm_only=false) {
        return getAllFullMatchingCoefficients(group, order, sm_only);
    }

    std::map<WCoef, complex_t> getAllFullRunCoefficients(WGroup group, QCDOrder order, bool sm_only=false) {
        std::vector<WCoef> ids = WCoefMapper::get_group(group);
        std::map<WCoef, complex_t> Cs;
        for (auto c : ids) {
            Cs.emplace(c, getFullRunCoefficient(group, c, order, sm_only));
        }

        return Cs;
    }

    std::map<WCoef, complex_t> getAFR(WGroup group, QCDOrder order, bool sm_only=false) {
        return getAllFullRunCoefficients(group, order, sm_only);
    }

    inline void build(WilsonConfig config) {
        std::map<std::string, std::shared_ptr<CoefficientGroup>> groups;
        auto model = ModelAPI().get();
       
        for (auto& gn : config.groups) {
            std::string gn_str = GroupMapper::str(gn);
            groups.emplace(gn_str, group_ptrs.at(gn_str));
            if (model == Model::THDM || model == Model::SUSY) {
                std::string gn_bsm_str = gn_str + "_" + ModelMapper::str(model);
                groups.emplace(gn_bsm_str, group_ptrs.at(gn_bsm_str));
            }
        }
        this->wm = CoefficientManager::Builder(ModelMapper::str(model), groups, config.matching_scale, config.hadronic_scale, OrderMapper::str(config.order));
    }
};

#endif // __WILSONINTERFACE_H__