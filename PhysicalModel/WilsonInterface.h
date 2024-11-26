#include "WilsonManager.h"
#include "Wilson_susyv2.h"
#include "Wilson_THDMv2.h"
#include "MartyWilson.h"
#include "General.h"

#include <map>



// enum class order {
//     LO, NLO, NNLO
// };

class WilsonInterface {

    CoefficientManager* wm;
    std::string model{"SM"};
    std::map<WilsonGroups, std::string> stringmapper = {{WilsonGroups::BCoefficients, "BCoefficients"}, {WilsonGroups::BPrimeCoefficients, "BPrimeCoefficients"}, {WilsonGroups::BScalarCoefficients, "BScalarCoefficients"}};
     std::map<std::string, std::shared_ptr<CoefficientGroup>> groupmapper 
     = {
        {"BCoefficients_SM", std::make_shared<BCoefficientGroup>()},
        {"BCoefficients_Prime_SM", std::make_shared<BPrimeCoefficientGroup>()},
        {"BCoefficients_Scalar_SM", std::make_shared<BScalarCoefficientGroup>()},
        {"BCoefficients_THDM", std::make_shared<BCoefficientGroup_THDM>()},
        {"BCoefficients_Prime_THDM", std::make_shared<BPrimeCoefficientGroup_THDM>()},
        {"BCoefficients_Scalar_THDM", std::make_shared<BScalarCoefficientGroup_THDM>()},
        {"BCoefficients_SUSY", std::make_shared<BCoefficientGroup_susy>()},
        {"BCoefficients_Prime_SUSY", std::make_shared<BPrimeCoefficientGroup_susy>()},
        {"BCoefficients_Scalar_SUSY", std::make_shared<BScalarCoefficientGroup_susy>()},
        {"BScalarCoefficients_Marty", std::make_shared<BCoefficientGroupMarty>()},
        {"BCoefficients_Prime_Marty", std::make_shared<BPrimeCoefficientGroupMarty>()},
        {"BCoefficients_Scalar_Marty", std::make_shared<BScalarCoefficientGroupMarty>()}};
    // };

    std::map<CoefficientOrder, std::string> ordermapping = {{CoefficientOrder::LO, "LO"}, {CoefficientOrder::NLO, "NLO"}, {CoefficientOrder::NNLO, "NNLO"}};
    std::string fake{"BCoefficients_SM"};
public:
    explicit WilsonInterface(const std::string& model) {
        this->wm = CoefficientManager::GetInstance(model);
    }

    void AddWilsonGroup(WilsonGroups groupname) {
        this->wm->registerCoefficientGroup(this->stringmapper[groupname], this->groupmapper[this->stringmapper[groupname]]);
        this->fake=this->stringmapper[groupname];
    }

    // ~WilsonInterface() {
    //     free(wm);
    // }

    void setQMatch(WilsonGroups groupName, double Q_match) {
        this->wm->setQMatch(this->stringmapper[groupName], Q_match);
    }

    void setParams(const std::string& block, int pdgCode, double value) {
        this->wm->setParams(fake, block, pdgCode, value);
    }

    void setMatchingCoefficient(WilsonGroups groupName, CoefficientOrder Order) {
        this->wm->setMatchingCoefficient(this->stringmapper[groupName], this->ordermapping[Order]);
    }

    void setGroupScale(WilsonGroups groupName, double Q) {
        this->wm->setGroupScale(this->stringmapper[groupName], Q);
    }

    void setRunCoefficient(WilsonGroups groupName, CoefficientOrder Order) {
        this->wm->setRunCoefficient(this->stringmapper[groupName], this->ordermapping[Order]);
    }

    void switchbasis(WilsonGroups groupName) {
        this->wm->switchbasis(this->stringmapper[groupName]);
    }

    double getAlphaS(WilsonGroups groupName) {
        return this->wm->getAlphaS(this->stringmapper[groupName]);
    }

    complex_t getMatchingCoefficient(WilsonGroups groupName, WilsonCoefficientList coeff, CoefficientOrder order) {
        return this->wm->getMatchingCoefficient(GroupMapper::str(groupName), WCoefMapper::str(coeff), OrderMapper::str(order));
    }

};