#if !defined(HYPERISO_WILSON_H)
#define HYPERISO_WILSON_H

#include "Include.h"
#include "Math.h"
#include "Utils.h"
#include "Wilson_parameters.h"
#include "ModelAPI.h"
#include "HasWilsonAPI.h"
#include "ParameterProxy.h"

class WilsonCoefficient {
public:
    WilsonCoefficient() = default;
    WilsonCoefficient(const std::string& name, const std::string& storage_block) : coeffName(name), storage_block(storage_block) {}

    // TODO : Implement initialization as dependent parameter from MARTY library or from lha
    void init(QCDOrder order);

    void set_name(std::string name) {this->coeffName = name;}
    void set_owned(bool owned);

    complex_t get_matching_value(std::string order) const; 
    std::string get_name() const {return this->coeffName;}
    QCDOrder get_max_order() const {return this->max_order;}
    LhaID id(QCDOrder order) const;

    bool operator==(const WilsonCoefficient& other) const;
    bool operator!=(const WilsonCoefficient& other) const { return !(*this == other); }

    virtual ~WilsonCoefficient() = default;

protected:
    virtual void LO_calculation() = 0;
    virtual void NLO_calculation() = 0;
    virtual void NNLO_calculation() = 0;

    std::string coeffName{};
    QCDOrder max_order;
    ContributionType type {ContributionType::SM};
    bool is_owned {false};
    bool from_lha {false};
    std::string storage_block;
};

// TODO : virtual LO, NLO, NNLO calculation need to be dealt properly

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
    os << "WilsonCoefficient " << coeff.get_name() << "has matching value: " << coeff.get_matching_value("LO") << " at LO" << std::endl;
    os << ", " << coeff.get_matching_value("NLO") << " at NLO, " << coeff.get_matching_value("NNLO") << "at NNLO" << std::endl;
    return os;
}



#endif //Wilsonv2