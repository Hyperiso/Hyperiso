#if !defined(HYPERISO_WILSON_H)
#define HYPERISO_WILSON_H

#include "Include.h"
#include "Math.h"
#include "Utils.h"
#include "Wilson_parameters.h"
#include "ModelAPI.h"
#include "HasWilsonAPI.h"
#include "ModelParamAdapter.h"

class WilsonCoefficient {
protected:
    WilsonCoefficient() {}

    void is_now_calculated(std::string order) {this->is_calculated[order] = true;}

    std::string CoeffName{};
    bool bsm {false};
    bool is_owned {false};

public:
    void set_name(std::string name) {this->CoeffName = name;}

    void set_WilsonCoeffMatching(const std::string& order, complex_t value) {
        std::cout << "Wow mais c'est encore plus génial ici" << std::endl; //TODO
    }

    void set_WilsonCoeffRun(const std::string& order, complex_t value) {
        std::cout << "Wow mais c'est vraiment encore plus génial ici" << std::endl; //TODO
    }

    // TODO : Store in parameters instead
    // void set_WilsonCoeffRun(std::string order, complex_t value) {this->CoefficientRunValue[order] = value;}
    // void set_WilsonCoeffMatching(std::string order ,complex_t value) {this->CoefficientMatchingValue[order] = value;}

    complex_t get_CoefficientMatchingValue(std::string order) const {
        auto base_id = WCoefMapper::flha_base(WCoefMapper::enum_elt(this->CoeffName));
        auto order_id = static_cast<long>(OrderMapper::enum_elt(order));
        LhaID code(base_id.first, base_id.second, order_id, static_cast<long>(this->bsm));
        ParameterProvider wilson_p = ParameterProvider(ParameterType::WILSON);
        return complex_t(wilson_p("B_MATCH", code), wilson_p("IMB_MATCH", code));
    }

    complex_t get_CoefficientRunValue(std::string order) const {return {0,0}; } //TODO

    complex_t get_CoefficientFullMatchingValue(std::string order) const {
        // double fact = QCDHelper::alpha_s(Q_match) / (4 * M_PI);

        // if (order == "LO") {
        //     return this->get_CoefficientMatchingValue("LO");
        // }
        // else if (order == "NLO") {
        //     return this->get_CoefficientMatchingValue("LO") + fact *  this->get_CoefficientMatchingValue("NLO");
        // }
        // else if (order == "NNLO") {
        //     return this->get_CoefficientMatchingValue("LO") + fact * this->get_CoefficientMatchingValue("NLO") + fact * fact * this->get_CoefficientMatchingValue("NNLO");
        // }
        // else {
        //     LOG_ERROR("ValueError", "Order request for wilson getfullmatching is not possible.");
        // }
    }


    complex_t get_CoefficientFullRunValue(std::string order) const {
        // double fact = QCDHelper::alpha_s(Q) / (4 * M_PI);

        // if (order == "LO") {
        //     return this->get_CoefficientRunValue("LO");
        // }
        // else if (order == "NLO") {
        //     return this->get_CoefficientRunValue("LO") + fact * this->get_CoefficientRunValue("NLO");
        // }
        // else if (order == "NNLO") {
        //     return this->get_CoefficientRunValue("LO") + fact *  this->get_CoefficientRunValue("NLO") + fact * fact * this->get_CoefficientRunValue("NNLO");
        // }
        // else {
        //     LOG_ERROR("ValueError", "Order request for wilson getfullrun is not possible.");
        // }
    }

    std::string get_name() const {return this->CoeffName;}

    virtual void LO_calculation() = 0;
    virtual void NLO_calculation() = 0;
    virtual void NNLO_calculation() = 0;

    bool is_it_calculated(std::string order) {return this->is_calculated[order];}

    void set_owned(bool owned);

    bool operator!=(const WilsonCoefficient& other) {
        return this->CoeffName != other.CoeffName;
    }
    WilsonCoefficient& operator+=(const WilsonCoefficient& other) {

        // if ((this->get_Q_match() == other.get_Q_match()) && (this->get_Q() == other.get_Q())) {
            //TODO
            // this->CoefficientMatchingValue["LO"] += other.get_CoefficientMatchingValue("LO");
            // this->CoefficientMatchingValue["NLO"] += other.get_CoefficientMatchingValue("NLO");
            // this->CoefficientMatchingValue["NNLO"] += other.get_CoefficientMatchingValue("NNLO");
            // this->CoefficientRunValue["LO"] += other.get_CoefficientRunValue("LO");
            // this->CoefficientRunValue["NLO"] += other.get_CoefficientRunValue("NLO");
            // this->CoefficientRunValue["NNLO"] += other.get_CoefficientRunValue("NNLO");
        // }
        // else {
        //     LOG_ERROR("ValueError", "Not the right scale"); 
        // }
        return *this;
    }

    virtual ~WilsonCoefficient() = default;

private:
    bool from_lha {false};
    std::map<std::string, bool> is_calculated{{"LO", false}, {"NLO", false}, {"NNLO", false}};
};

//TODO : virtual LO,NLO,NNLO calculation need to be dealed properly

class C_Blnu_A : public WilsonCoefficient {
public:
    C_Blnu_A() : WilsonCoefficient() {this->set_name("C_Blnu_A");}

    void LO_calculation() { }; 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_Blnu_P : public WilsonCoefficient {
public:
    C_Blnu_P() : WilsonCoefficient() {this->set_name("C_Blnu_P");}

    void LO_calculation() {}
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_V1 : public WilsonCoefficient {
public:
    C_V1() : WilsonCoefficient() {this->set_name("C_V1");}

    void LO_calculation() {
         //return double_to_complex_save("LO", 1.);  TODO
        } 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_V2 : public WilsonCoefficient {
public:
    C_V2() : WilsonCoefficient() {this->set_name("C_V2");}

    void LO_calculation() {}
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_S1 : public WilsonCoefficient {
public:
    C_S1() : WilsonCoefficient() {this->set_name("C_S1");}

    void LO_calculation() {}
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_S2 : public WilsonCoefficient {
public:
    C_S2() : WilsonCoefficient() {this->set_name("C_S2");}

    void LO_calculation() {}
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

class C_T : public WilsonCoefficient {
public:
    C_T() : WilsonCoefficient() {this->set_name("C_T");}

    void LO_calculation() {}
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

};

inline std::ostream& operator<<(std::ostream& os, WilsonCoefficient& coeff) {
    os << "WilsonCoefficient " << coeff.get_name() << "has matching value: " << coeff.get_CoefficientMatchingValue("LO") << " at LO" << std::endl;
    os<< ", " << coeff.get_CoefficientMatchingValue("NLO") << " at NLO, " << coeff.get_CoefficientMatchingValue("NNLO") << "at NNLO" << std::endl;
    os << "and " << "has run value: " << coeff.get_CoefficientRunValue("LO") << " at LO" << std::endl;
    os<< ", " << coeff.get_CoefficientRunValue("NLO") << " at NLO, " << coeff.get_CoefficientRunValue("NNLO") << "at NNLO" << std::endl;
    return os;
}



#endif //Wilsonv2