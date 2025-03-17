#ifndef __WILSONINTERFACE_H__
#define __WILSONINTERFACE_H__

#include "MemoryManager.h"
#include "WilsonManager.h"
#ifdef BUILD_WITH_SOFTSUSY
#include "Wilson_SUSY.h"
#endif
#ifdef BUILD_WITH_2HDMC
#include "Wilson_THDM.h"
#endif
#include "MartyWilson.h"
#include "General.h"

#include <map>

class WilsonInterface {
private:
    std::shared_ptr<CoefficientManager> wm;

    std::map<std::string, std::shared_ptr<CoefficientGroup>> group_ptrs = {
        {"BCoefficients", std::make_shared<BCoefficientGroup>()},
        {"BPrimeCoefficients", std::make_shared<BPrimeCoefficientGroup>()},
        {"BScalarCoefficients", std::make_shared<BScalarCoefficientGroup>()},
        {"BlnuCoefficients", std::make_shared<BlnuCoefficientGroup>()},
        {"BclnuCoefficients", std::make_shared<BclnuCoefficientGroup>()},
        #ifdef BUILD_WITH_2HDMC
        {"BCoefficients_THDM", std::make_shared<BCoefficientGroup_THDM>()},
        {"BPrimeCoefficients_THDM", std::make_shared<BPrimeCoefficientGroup_THDM>()},
        {"BScalarCoefficients_THDM", std::make_shared<BScalarCoefficientGroup_THDM>()},
        {"BlnuCoefficients_THDM", std::make_shared<BlnuCoefficientGroup_THDM>()},
        {"BclnuCoefficients_THDM", std::make_shared<BclnuCoefficientGroup_THDM>()},
        #endif
        #ifdef BUILD_WITH_SOFTSUSY
        {"BCoefficients_SUSY", std::make_shared<BCoefficientGroup_susy>()},
        {"BPrimeCoefficients_SUSY", std::make_shared<BPrimeCoefficientGroup_susy>()},
        {"BScalarCoefficients_SUSY", std::make_shared<BScalarCoefficientGroup_susy>()},
        {"BlnuCoefficients_SUSY", std::make_shared<BlnuCoefficientGroup_SUSY>()},
        {"BclnuCoefficients_SUSY", std::make_shared<BclnuCoefficientGroup_SUSY>()},
        #endif
    };

    QCDOrder ensure_mty_compat(QCDOrder order) {
        if (MemoryManager::GetInstance()->useMarty() && !(order == QCDOrder::LO)) {
            LOG_WARN("Using MARTY defaults all calculations to LO in QCD.");
            return QCDOrder::LO;
        }
        return order;
    }

public:
    explicit WilsonInterface() {
        this->wm = CoefficientManager::GetInstance();
    }

    void addWilsonGroup(WGroup group_name) {
        this->wm->registerCoefficientGroup(GroupMapper::str(group_name), this->group_ptrs.at(GroupMapper::str(group_name)));
    }

    void setQMatch(WGroup group_name, double Q_match) {
        this->wm->setQMatch(GroupMapper::str(group_name), Q_match);
    }

    void setParams(const std::string& block, int pdgCode, double value) {
        this->wm->setParams("BCoefficients", block, pdgCode, value);
    }

    void setMatchingCoefficient(WGroup group_name, QCDOrder order) {
        this->wm->setMatchingCoefficient(GroupMapper::str(group_name), OrderMapper::str(ensure_mty_compat(order)));
    }

    void setGroupScale(WGroup group_name, double Q) {
        this->wm->setGroupScale(GroupMapper::str(group_name), Q);
    }

    void setRunCoefficient(WGroup group_name, QCDOrder order) {
        this->wm->setRunCoefficient(GroupMapper::str(group_name), OrderMapper::str(ensure_mty_compat(order)));
    }

    void switchbasis(WGroup group_name) {
        this->wm->switchbasis(GroupMapper::str(group_name));
    }

    double getAlphaS(WGroup group_name) {
        return this->wm->getAlphaS(GroupMapper::str(group_name));
    }

    complex_t getMatchingCoefficient(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return this->wm->getMatchingCoefficient(GroupMapper::str(group), WCoefMapper::str(coeff), OrderMapper::str(ensure_mty_compat(order)), sm_only);
    }

    complex_t getM(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return getMatchingCoefficient(group, coeff, order, sm_only);
    }

    complex_t getFullMatchingCoefficient(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return this->wm->getFullMatchingCoefficient(GroupMapper::str(group), WCoefMapper::str(coeff), OrderMapper::str(ensure_mty_compat(order)), sm_only);
    }

    complex_t getFM(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return getFullMatchingCoefficient(group, coeff, order, sm_only);
    }

    complex_t getRunCoefficient(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return this->wm->getRunCoefficient(GroupMapper::str(group), WCoefMapper::str(coeff), OrderMapper::str(ensure_mty_compat(order)), sm_only);
    }

    complex_t getR(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return getRunCoefficient(group, coeff, order, sm_only);
    }

    complex_t getFullRunCoefficient(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) {
        return this->wm->getFullRunCoefficient(GroupMapper::str(group), WCoefMapper::str(coeff), OrderMapper::str(ensure_mty_compat(order)), sm_only);
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

    inline void build(std::vector<WGroup> group_names, double Q_match, double Q, QCDOrder order) {
        std::map<std::string, std::shared_ptr<CoefficientGroup>> groups;
        auto model = MemoryManager::GetInstance()->getModel();
        for (auto& gn : group_names) {
            std::string gn_str = GroupMapper::str(gn);
            groups.emplace(gn_str, group_ptrs.at(gn_str));
            if (model == Model::THDM || model == Model::SUSY) {
                std::string gn_bsm_str = gn_str + "_" + ModelMapper::str(model);
                groups.emplace(gn_bsm_str, group_ptrs.at(gn_bsm_str));
            }
        }
        this->wm = CoefficientManager::Builder(ModelMapper::str(model), groups, Q_match, Q, OrderMapper::str(order));
    }
};

#endif // __WILSONINTERFACE_H__