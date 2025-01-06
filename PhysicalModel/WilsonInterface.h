#ifndef __WILSONINTERFACE_H__
#define __WILSONINTERFACE_H__

#include "WilsonManager.h"
#ifdef BUILD_WITH_SOFTSUSY
#include "Wilson_susyv2.h"
#endif
#ifdef BUILD_WITH_2HDMC
#include "Wilson_THDMv2.h"
#endif
#include "MartyWilson.h"
#include "General.h"

#include <map>

class WilsonInterface {

    std::shared_ptr<CoefficientManager> wm;
    std::string model{"SM"};
     std::map<std::string, std::shared_ptr<CoefficientGroup>> groupmapper 
     = {
        {"BCoefficients_SM", std::make_shared<BCoefficientGroup>()},
        {"BCoefficients_Prime_SM", std::make_shared<BPrimeCoefficientGroup>()},
        {"BCoefficients_Scalar_SM", std::make_shared<BScalarCoefficientGroup>()},
        #ifdef BUILD_WITH_2HDMC
        {"BCoefficients_THDM", std::make_shared<BCoefficientGroup_THDM>()},
        {"BCoefficients_Prime_THDM", std::make_shared<BPrimeCoefficientGroup_THDM>()},
        {"BCoefficients_Scalar_THDM", std::make_shared<BScalarCoefficientGroup_THDM>()},
        #endif
        #ifdef BUILD_WITH_SOFTSUSY
        {"BCoefficients_SUSY", std::make_shared<BCoefficientGroup_susy>()},
        {"BCoefficients_Prime_SUSY", std::make_shared<BPrimeCoefficientGroup_susy>()},
        {"BCoefficients_Scalar_SUSY", std::make_shared<BScalarCoefficientGroup_susy>()},
        #endif
        {"BScalarCoefficients_Marty", std::make_shared<BCoefficientGroupMarty>()},
        {"BCoefficients_Prime_Marty", std::make_shared<BPrimeCoefficientGroupMarty>()},
        {"BCoefficients_Scalar_Marty", std::make_shared<BScalarCoefficientGroupMarty>()}
        };
    // };

    std::map<QCDOrder, std::string> ordermapping = {{QCDOrder::LO, "LO"}, {QCDOrder::NLO, "NLO"}, {QCDOrder::NNLO, "NNLO"}};
    std::string fake{"BCoefficients_SM"};
public:
    explicit WilsonInterface(const std::string& model) {
        this->wm = CoefficientManager::GetInstance(model);
    }

    void AddWilsonGroup(WilsonGroups groupname) {
        this->wm->registerCoefficientGroup(GroupMapper::str(groupname), this->groupmapper[GroupMapper::str(groupname)]);
        this->fake=GroupMapper::str(groupname);
    }


    void setQMatch(WilsonGroups groupName, double Q_match) {
        this->wm->setQMatch(GroupMapper::str(groupName), Q_match);
    }

    void setParams(const std::string& block, int pdgCode, double value) {
        this->wm->setParams(fake, block, pdgCode, value);
    }

    void setMatchingCoefficient(WilsonGroups groupName, QCDOrder Order) {
        if (MemoryManager::GetInstance()->useMarty() && !(Order == QCDOrder::LO)) {
            LOG_WARN("Careful, cannot go above LO using Marty for now.");
            Order = QCDOrder::LO;
        }
        this->wm->setMatchingCoefficient(GroupMapper::str(groupName), this->ordermapping[Order]);
    }

    void setGroupScale(WilsonGroups groupName, double Q) {
        this->wm->setGroupScale(GroupMapper::str(groupName), Q);
    }

    void setRunCoefficient(WilsonGroups groupName, QCDOrder Order) {
        if (MemoryManager::GetInstance()->useMarty() && !(Order == QCDOrder::LO)) {
            LOG_WARN("Careful, cannot go above LO using Marty for now.");
            Order = QCDOrder::LO;
        }
        this->wm->setRunCoefficient(GroupMapper::str(groupName), this->ordermapping[Order]);
    }

    void switchbasis(WilsonGroups groupName) {
        this->wm->switchbasis(GroupMapper::str(groupName));
    }

    double getAlphaS(WilsonGroups groupName) {
        return this->wm->getAlphaS(GroupMapper::str(groupName));
    }

    complex_t getMatchingCoefficient(WilsonGroups groupName, BWilsonCoefficients coeff, QCDOrder order) {
        if (MemoryManager::GetInstance()->useMarty() && !(order == QCDOrder::LO)) {
            LOG_WARN("We only have LO contribution while using Marty, please set order to LO.");
            order = QCDOrder::LO;
        }
        return this->wm->getMatchingCoefficient(GroupMapper::str(groupName), WCoefMapper::str(coeff), OrderMapper::str(order));
    }

    inline void Builder(std::string model, std::vector<WilsonGroups> grp_list, double Q_match, double Q, QCDOrder order) {
        std::map<std::string, std::shared_ptr<CoefficientGroup>> groups;
        for (auto& wil : grp_list) {
            groups[GroupMapper::str(wil)] = groupmapper[GroupMapper::str(wil)];
        }
        this->wm = CoefficientManager::Builder(model, groups, Q_match, Q, OrderMapper::str(order));
    }

};

#endif // __WILSONINTERFACE_H__