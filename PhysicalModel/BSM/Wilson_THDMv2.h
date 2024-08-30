#if !defined(HYPERISO_WILSON_THDM_H)
#define HYPERISO_WILSON_THDM_H
#include "Wilsonv2.h"
#include "thdm_parameters.h"
#include "Math_THDM.h"

class WilsonCoefficient_THDM {
protected:
    WilsonCoefficient_THDM(double Q_match) {thdm_params->set_params(Q_match);}
    WilsonCoefficient_THDM() {thdm_params->set_params(81.);}
    void set_mod_parameters(Parameters* new_mod) {this->mod = new_mod;};
    

    Parameters* mod = Parameters::GetInstance(2);

    thdm_parameters *thdm_params = thdm_parameters::GetInstance();
};
class C1_THDM : public C1, public WilsonCoefficient_THDM {
public:
    C1_THDM(double Q_match) : C1(Q_match) {this->set_name("C1_THDM");}
    C1_THDM(double Q_match, int gen) : C1(Q_match) {}
    C1_THDM() : C1() {this->set_name("C1_THDM");

    }

    std::complex<double> LO_calculation() { this->set_CoefficientMatchingValue("LO", {0.,0.}); std::cout << "coefficient C1 " << this->get_CoefficientMatchingValue("LO") << std::endl; return {0,0};} 
    std::complex<double> NLO_calculation() { this->set_CoefficientMatchingValue("NLO", {0.,0.}); std::cout << "coefficient C1 NLO " << this->get_CoefficientMatchingValue("NLO") << std::endl;return {0,0};} 
    std::complex<double> NNLO_calculation() { this->set_CoefficientMatchingValue("NNLO", {0.,0.}); std::cout << "coefficient C1 NNLO " << this->get_CoefficientMatchingValue("NNLO") << std::endl;return {0,0};} 

};

class C2_THDM : public C2, public WilsonCoefficient_THDM {
public:
    C2_THDM(double Q_match) : C2(Q_match), WilsonCoefficient_THDM(Q_match){}
    C2_THDM(double Q_match, int gen) : C2(Q_match), WilsonCoefficient_THDM(Q_match) {}
    
    C2_THDM() : C2() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation() {this->set_CoefficientMatchingValue("NLO", {0.,0.}); return {0,0};} 
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class C3_THDM : public C3, public WilsonCoefficient_THDM {
public:
    C3_THDM(double Q_match) : C3(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C3_THDM(double Q_match, int gen) : C3(Q_match), WilsonCoefficient_THDM(Q_match) {}
    
    C3_THDM() : C3() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation()  {this->set_CoefficientMatchingValue("NLO", {0.,0.}); return {0,0};} 
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class C4_THDM : public C4, public WilsonCoefficient_THDM {
public:
    C4_THDM(double Q_match) : C4(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C4_THDM(double Q_match, int gen) : C4(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C4_THDM() : C4() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class C5_THDM : public C5, public WilsonCoefficient_THDM {
public:
    C5_THDM(double Q_match) : C5(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C5_THDM(double Q_match, int gen) : C5(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C5_THDM() : C5() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation() {this->set_CoefficientMatchingValue("NLO", {0.,0.}); return {0,0};} 
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class C6_THDM : public C6, public WilsonCoefficient_THDM {
public:
    C6_THDM(double Q_match) : C6(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C6_THDM(double Q_match, int gen) : C6(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C6_THDM() : C6() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation() {this->set_CoefficientMatchingValue("NLO", {0.,0.}); return {0,0};} 
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class C7_THDM : public C7, public WilsonCoefficient_THDM {
public:
    C7_THDM(double Q_match) : C7(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C7_THDM(double Q_match, int gen) : C7(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C7_THDM() : C7() {

    }

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class C8_THDM : public C8, public WilsonCoefficient_THDM {
public:
    C8_THDM(double Q_match) : C8(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C8_THDM(double Q_match, int gen) : C8(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C8_THDM() : C8() {

    }

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class C9_THDM : public C9, public WilsonCoefficient_THDM {
public:
    C9_THDM(double Q_match) : C9(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C9_THDM(double Q_match, int gen) : C9(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C9_THDM() : C9() {

    }

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class C10_THDM : public C10, public WilsonCoefficient_THDM {
public:
    C10_THDM(double Q_match) : C10(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C10_THDM(double Q_match, int gen) : C10(Q_match), WilsonCoefficient_THDM(Q_match) {}
    C10_THDM() : C10() {

    }

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CQ1_THDM : public CQ1, public WilsonCoefficient_THDM {
public:
    CQ1_THDM(double Q_match) : CQ1(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CQ1_THDM(double Q_match, int gen) : CQ1(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CQ1_THDM() : CQ1() {

    }

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation() {this->set_CoefficientMatchingValue("NLO", {0.,0.}); return {0,0};}
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CQ2_THDM : public CQ2, public WilsonCoefficient_THDM {
public:
    CQ2_THDM(double Q_match) : CQ2(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CQ2_THDM(double Q_match, int gen) : CQ2(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CQ2_THDM() : CQ2() {

    }

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation() {this->set_CoefficientMatchingValue("NLO", {0.,0.}); return {0,0};}
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CP1_THDM : public CP1, public WilsonCoefficient_THDM {
public:
    CP1_THDM(double Q_match) : CP1(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP1_THDM(double Q_match, int gen) : CP1(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP1_THDM() : CP1() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CP2_THDM : public CP2, public WilsonCoefficient_THDM {
public:
    CP2_THDM(double Q_match) : CP2(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP2_THDM(double Q_match, int gen) : CP2(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP2_THDM() : CP2() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CP3_THDM : public CP3, public WilsonCoefficient_THDM {
public:
    CP3_THDM(double Q_match) : CP3(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP3_THDM(double Q_match, int gen) : CP3(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP3_THDM() : CP3() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CP4_THDM : public CP4, public WilsonCoefficient_THDM {
public:
    CP4_THDM(double Q_match) : CP4(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP4_THDM(double Q_match, int gen) : CP4(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP4_THDM() : CP4() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CP5_THDM : public CP5, public WilsonCoefficient_THDM {
public:
    CP5_THDM(double Q_match) : CP5(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP5_THDM(double Q_match, int gen) : CP5(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP5_THDM() : CP5() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CP6_THDM : public CP6, public WilsonCoefficient_THDM {
public:
    CP6_THDM(double Q_match) : CP6(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP6_THDM(double Q_match, int gen) : CP6(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP6_THDM() : CP6() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CP7_THDM : public CP7, public WilsonCoefficient_THDM {
public:
    CP7_THDM(double Q_match) : CP7(Q_match) , WilsonCoefficient_THDM(Q_match){}
    CP7_THDM(double Q_match, int gen) : CP7(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP7_THDM() : CP7() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CP8_THDM : public CP8, public WilsonCoefficient_THDM {
public:
    CP8_THDM(double Q_match) : CP8(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP8_THDM(double Q_match, int gen) : CP8(Q_match) , WilsonCoefficient_THDM(Q_match) {}
    CP8_THDM() : CP8() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CP9_THDM : public CP9, public WilsonCoefficient_THDM {
public:
    CP9_THDM(double Q_match) : CP9(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP9_THDM(double Q_match, int gen) : CP9(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP9_THDM() : CP9() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CP10_THDM : public CP10, public WilsonCoefficient_THDM {
public:
    CP10_THDM(double Q_match) : CP10(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP10_THDM(double Q_match, int gen) : CP10(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CP10_THDM() : CP10() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NLO", {0.,0.}); return {0,0};} 

};

class CPQ1_THDM : public CPQ1, public WilsonCoefficient_THDM {
public:
    CPQ1_THDM(double Q_match) : CPQ1(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CPQ1_THDM(double Q_match, int gen) : CPQ1(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CPQ1_THDM() : CPQ1() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class CPQ2_THDM : public CPQ2, public WilsonCoefficient_THDM {
public:
    CPQ2_THDM(double Q_match) : CPQ2(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CPQ2_THDM(double Q_match, int gen) : CPQ2(Q_match), WilsonCoefficient_THDM(Q_match) {}
    CPQ2_THDM() : CPQ2() {

    }

    std::complex<double> LO_calculation() {this->set_CoefficientMatchingValue("LO", {0.,0.}); return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {this->set_CoefficientMatchingValue("NNLO", {0.,0.}); return {0,0};} 

};

class BCoefficientGroup_THDM : public BCoefficientGroup {

public:
    BCoefficientGroup_THDM() {this->clear();
        this->insert(std::make_pair("C1", std::make_unique<C1_THDM>())); this->insert(std::make_pair("C2", std::make_unique<C2_THDM>())); this->insert(std::make_pair("C3", std::make_unique<C3_THDM>()));
        this->insert(std::make_pair("C4", std::make_unique<C4_THDM>()));  this->insert(std::make_pair("C5", std::make_unique<C5_THDM>())); this->insert(std::make_pair("C6", std::make_unique<C6_THDM>())); 
        this->insert(std::make_pair("C7", std::make_unique<C7_THDM>()));  this->insert(std::make_pair("C8", std::make_unique<C8_THDM>()));  this->insert(std::make_pair("C9", std::make_unique<C9_THDM>())); 
        this->insert(std::make_pair("C10", std::make_unique<C10_THDM>())); 
    }
    BCoefficientGroup_THDM(double Q_match) {this->clear();
        this->insert(std::make_pair("C1", std::make_unique<C1_THDM>(Q_match))); this->insert(std::make_pair("C2", std::make_unique<C2_THDM>(Q_match))); this->insert(std::make_pair("C3", std::make_unique<C3_THDM>(Q_match)));
        this->insert(std::make_pair("C4", std::make_unique<C4_THDM>(Q_match)));  this->insert(std::make_pair("C5", std::make_unique<C5_THDM>(Q_match))); this->insert(std::make_pair("C6", std::make_unique<C6_THDM>(Q_match))); 
        this->insert(std::make_pair("C7", std::make_unique<C7_THDM>(Q_match)));  this->insert(std::make_pair("C8", std::make_unique<C8_THDM>(Q_match)));  this->insert(std::make_pair("C9", std::make_unique<C9_THDM>(Q_match))); 
        this->insert(std::make_pair("C10", std::make_unique<C10_THDM>(Q_match)));
    }
    BCoefficientGroup_THDM(double Q_match, int gen) {
        this->insert(std::make_pair("C1", std::make_unique<C1_THDM>(Q_match, gen))); this->insert(std::make_pair("C2", std::make_unique<C2_THDM>(Q_match, gen))); this->insert(std::make_pair("C3", std::make_unique<C3_THDM>(Q_match, gen)));
        this->insert(std::make_pair("C4", std::make_unique<C4_THDM>(Q_match, gen)));  this->insert(std::make_pair("C5", std::make_unique<C5_THDM>(Q_match, gen))); this->insert(std::make_pair("C6", std::make_unique<C6_THDM>(Q_match, gen))); 
        this->insert(std::make_pair("C7", std::make_unique<C7_THDM>(Q_match, gen)));  this->insert(std::make_pair("C8", std::make_unique<C8_THDM>(Q_match, gen)));  this->insert(std::make_pair("C9", std::make_unique<C9_THDM>(Q_match, gen))); 
        this->insert(std::make_pair("C10", std::make_unique<C10_THDM>(Q_match, gen)));
    }

};

class BPrimeCoefficientGroup_THDM : public BPrimeCoefficientGroup {
public:
    BPrimeCoefficientGroup_THDM() {
        this->insert(std::make_pair("CP1", std::make_unique<CP1_THDM>())); this->insert(std::make_pair("CP2", std::make_unique<CP2_THDM>())); this->insert(std::make_pair("CP3", std::make_unique<CP3_THDM>()));
        this->insert(std::make_pair("CP4", std::make_unique<CP4_THDM>()));  this->insert(std::make_pair("CP5", std::make_unique<CP5_THDM>())); this->insert(std::make_pair("CP6", std::make_unique<CP6_THDM>())); 
        this->insert(std::make_pair("CP7", std::make_unique<CP7_THDM>()));  this->insert(std::make_pair("CP8", std::make_unique<CP8_THDM>()));  this->insert(std::make_pair("CP9", std::make_unique<CP9_THDM>())); 
        this->insert(std::make_pair("CP10", std::make_unique<CP10_THDM>())); this->insert(std::make_pair("CPQ1", std::make_unique<CPQ1_THDM>())); this->insert(std::make_pair("CPQ2", std::make_unique<CPQ2_THDM>())); 
    }
    BPrimeCoefficientGroup_THDM(double Q_match) {
        this->insert(std::make_pair("CP1", std::make_unique<CP1_THDM>(Q_match))); this->insert(std::make_pair("CP2", std::make_unique<CP2_THDM>(Q_match))); this->insert(std::make_pair("CP3", std::make_unique<CP3_THDM>(Q_match)));
        this->insert(std::make_pair("CP4", std::make_unique<CP4_THDM>(Q_match)));  this->insert(std::make_pair("CP5", std::make_unique<CP5_THDM>(Q_match))); this->insert(std::make_pair("CP6", std::make_unique<CP6_THDM>(Q_match))); 
        this->insert(std::make_pair("CP7", std::make_unique<CP7_THDM>(Q_match)));  this->insert(std::make_pair("CP8", std::make_unique<CP8_THDM>(Q_match)));  this->insert(std::make_pair("CP9", std::make_unique<CP9_THDM>(Q_match))); 
        this->insert(std::make_pair("CP10", std::make_unique<CP10_THDM>(Q_match))); this->insert(std::make_pair("CPQ1", std::make_unique<CPQ1_THDM>(Q_match))); this->insert(std::make_pair("CPQ2", std::make_unique<CPQ2_THDM>(Q_match)));
    }

    void set_base_1();
    void set_base_2();

};

class BScalarCoefficientGroup_THDM : public BScalarCoefficientGroup {
public:
    BScalarCoefficientGroup_THDM() {
        this->insert(std::make_pair("CQ1", std::make_unique<CQ1_THDM>())); this->insert(std::make_pair("CQ2", std::make_unique<CQ2_THDM>()));
    }
    BScalarCoefficientGroup_THDM(double Q_match) {
        this->insert(std::make_pair("CQ1", std::make_unique<CQ1_THDM>(Q_match))); this->insert(std::make_pair("CQ2", std::make_unique<CQ2_THDM>(Q_match)));
    }

    void set_base_1();
    void set_base_2();

};


#endif