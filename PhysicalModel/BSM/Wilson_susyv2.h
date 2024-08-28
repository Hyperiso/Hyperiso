#if !defined(HYPERISO_WILSON_SUSY_H)
#define HYPERISO_WILSON_SUSY_H
#include "Wilsonv2.h"
#include "susy_parameters.h"
#include "Math_THDM.h"

class WilsonCoefficient_susy {
protected:
    void set_mod_parameters(Parameters* new_mod) {this->mod = new_mod;};

    WilsonCoefficient_susy (double Q_match) {this->Q_match = Q_match; sus_param = susy_parameters::GetInstance(this->Q_match);
    calculateContribution = [&](auto hFunc, const Array3D_3x7x4& X, const Array3D_3x7x4& X2, int ie, int ae, bool isChargeps) -> double {
            double ratio = std::pow((*Parameters::GetInstance(0))("MASS", 24) / (*sus_param).Mch[ie], 2);
            double msqOverMchSquared = std::pow((*sus_param).MsqU[ae] / (*sus_param).Mch[ie], 2.0);
            double factor = isChargeps ? (-(*sus_param).epsilonb / (1.0 + (*sus_param).epsilonb * (*susy)("HMIX",2)) * (*susy)("HMIX",2)) : 1.0;
            return ratio * (
                X[ie][ae][1] * X2[ie][ae][2] * hFunc(msqOverMchSquared)) * (*sus_param).kappaFactor * factor;
        };
    }
    
    Parameters* mod = Parameters::GetInstance(2);
    int gen{3};
    double Q_match;
    
    Parameters* susy = Parameters::GetInstance(1);
	std::function<double(std::function<double(double)>, const Array3D_3x7x4&, const Array3D_3x7x4&, int, int, bool)> calculateContribution;

	susy_parameters* sus_param;
};
class C1_susy : public C1, public WilsonCoefficient_susy {
public:
    C1_susy() : C1(), WilsonCoefficient_susy(81.) {}
    C1_susy(double Q_match) : C1(Q_match), WilsonCoefficient_susy(Q_match) {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation()  {return {0,0};} 
    std::complex<double> NNLO_calculation();

};

class C2_susy : public C2, public WilsonCoefficient_susy {
public:
    C2_susy() : C2(), WilsonCoefficient_susy(81.) {}
    C2_susy(double Q_match) : C2(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation()  {return {0,0};} 
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C3_susy : public C3, public WilsonCoefficient_susy {
public:
    C3_susy() : C3(), WilsonCoefficient_susy(81.) {}
    C3_susy(double Q_match) : C3(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation() {return {0,0};} 
    std::complex<double> NNLO_calculation();

};

class C4_susy : public C4, public WilsonCoefficient_susy {
public:
    C4_susy() : C4(), WilsonCoefficient_susy(81.) {}
    C4_susy(double Q_match) : C4(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation(); 

};

class C5_susy : public C5, public WilsonCoefficient_susy {
public:
    C5_susy() : C5(), WilsonCoefficient_susy(81.) {}
    C5_susy(double Q_match) : C5(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation() {return {0,0};} 
    std::complex<double> NNLO_calculation();

};

class C6_susy : public C6, public WilsonCoefficient_susy {
public:
    C6_susy() : C6(), WilsonCoefficient_susy(81.) {}
    C6_susy(double Q_match) : C6(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation() {return {0,0};} 
    std::complex<double> NNLO_calculation();

};

class C7_susy : public C7, public WilsonCoefficient_susy {
public:
    C7_susy() : C7(), WilsonCoefficient_susy(81.) {}
    C7_susy(double Q_match) : C7(Q_match), WilsonCoefficient_susy(Q_match)  {}

    std::complex<double> LO_calculation() override;
    std::complex<double> NLO_calculation() override;
    std::complex<double> NNLO_calculation() override {return {0,0};} 

};

class C8_susy : public C8, public WilsonCoefficient_susy {
public:
    C8_susy() : C8(), WilsonCoefficient_susy(81.) {}
    C8_susy(double Q_match) : C8(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation()override;
    std::complex<double> NLO_calculation() override;
    std::complex<double> NNLO_calculation() override {return {0,0};} 

};

class C9_susy : public C9, public WilsonCoefficient_susy {
public:
    C9_susy() : C9(), WilsonCoefficient_susy(81.) {}
    C9_susy(double Q_match) : C9(Q_match), WilsonCoefficient_susy(Q_match)  {}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C10_susy : public C10, public WilsonCoefficient_susy {
public:
    C10_susy() : C10(), WilsonCoefficient_susy(81.) {}
    C10_susy(double Q_match) : C10(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CQ1_susy : public CQ1, public WilsonCoefficient_susy {
public:
    CQ1_susy() : CQ1(), WilsonCoefficient_susy(81.) {}
    CQ1_susy(double Q_match) : CQ1(Q_match), WilsonCoefficient_susy(Q_match) {}
    CQ1_susy(double Q_match, int gen) : CQ1(Q_match), WilsonCoefficient_susy(Q_match) {}


    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CQ2_susy : public CQ2, public WilsonCoefficient_susy {
public:
    CQ2_susy() : CQ2(), WilsonCoefficient_susy(81.) {}
    CQ2_susy(double Q_match) : CQ2(Q_match), WilsonCoefficient_susy(Q_match) {}
    CQ2_susy(double Q_match, int gen) : CQ2(Q_match), WilsonCoefficient_susy(Q_match) {}


    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP1_susy : public CP1, public WilsonCoefficient_susy {
public:
    CP1_susy() : CP1(), WilsonCoefficient_susy(81.) {}
    CP1_susy(double Q_match) : CP1(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP1_susy(double Q_match, int gen) : CP1(Q_match), WilsonCoefficient_susy(Q_match) {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP2_susy : public CP2, public WilsonCoefficient_susy {
public:
    CP2_susy() : CP2(), WilsonCoefficient_susy(81.) {}
    CP2_susy(double Q_match) : CP2(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP2_susy(double Q_match, int gen) : CP2(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP3_susy : public CP3, public WilsonCoefficient_susy {
public:
    CP3_susy() : CP3(), WilsonCoefficient_susy(81.) {}
    CP3_susy(double Q_match) : CP3(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP3_susy(double Q_match, int gen) : CP3(Q_match), WilsonCoefficient_susy(Q_match) {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP4_susy : public CP4, public WilsonCoefficient_susy {
public:
    CP4_susy() : CP4(), WilsonCoefficient_susy(81.) {}
    CP4_susy(double Q_match) : CP4(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP4_susy(double Q_match, int gen) : CP4(Q_match), WilsonCoefficient_susy(Q_match) {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP5_susy : public CP5, public WilsonCoefficient_susy {
public:
    CP5_susy() : CP5(), WilsonCoefficient_susy(81.) {}
    CP5_susy(double Q_match) : CP5(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP5_susy(double Q_match, int gen) : CP5(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP6_susy : public CP6, public WilsonCoefficient_susy {
public:
    CP6_susy() : CP6(), WilsonCoefficient_susy(81.) {}
    CP6_susy(double Q_match) : CP6(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP6_susy(double Q_match, int gen) : CP6(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP7_susy : public CP7, public WilsonCoefficient_susy {
public:
    CP7_susy() : CP7(), WilsonCoefficient_susy(81.) {}
    CP7_susy(double Q_match) : CP7(Q_match) , WilsonCoefficient_susy(Q_match){}
    CP7_susy(double Q_match, int gen) : CP7(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP8_susy : public CP8, public WilsonCoefficient_susy {
public:
    CP8_susy() : CP8(), WilsonCoefficient_susy(81.) {}
    CP8_susy(double Q_match) : CP8(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP8_susy(double Q_match, int gen) : CP8(Q_match) , WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation(); 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP9_susy : public CP9, public WilsonCoefficient_susy {
public:
    CP9_susy() : CP9(), WilsonCoefficient_susy(81.) {}
    CP9_susy(double Q_match) : CP9(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP9_susy(double Q_match, int gen) : CP9(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP10_susy : public CP10, public WilsonCoefficient_susy {
public:
    CP10_susy() : CP10(), WilsonCoefficient_susy(81.) {}
    CP10_susy(double Q_match) : CP10(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP10_susy(double Q_match, int gen) : CP10(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CPQ1_susy : public CPQ1, public WilsonCoefficient_susy {
public:
    CPQ1_susy() : CPQ1(), WilsonCoefficient_susy(81.) {}
    CPQ1_susy(double Q_match) : CPQ1(Q_match), WilsonCoefficient_susy(Q_match) {}
    CPQ1_susy(double Q_match, int gen) : CPQ1(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CPQ2_susy : public CPQ2, public WilsonCoefficient_susy {
public:
    CPQ2_susy() : CPQ2(), WilsonCoefficient_susy(81.) {}
    CPQ2_susy(double Q_match) : CPQ2(Q_match), WilsonCoefficient_susy(Q_match) {}
    CPQ2_susy(double Q_match, int gen) : CPQ2(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class BCoefficientGroup_susy : public BCoefficientGroup {

public:
    BCoefficientGroup_susy() {
        this->insert(std::make_pair("C1", std::make_unique<C1_susy>())); this->insert(std::make_pair("C2", std::make_unique<C2_susy>())); this->insert(std::make_pair("C3", std::make_unique<C3_susy>()));
        this->insert(std::make_pair("C4", std::make_unique<C4_susy>()));  this->insert(std::make_pair("C5", std::make_unique<C5_susy>())); this->insert(std::make_pair("C6", std::make_unique<C6_susy>())); 
        this->insert(std::make_pair("C7", std::make_unique<C7_susy>()));  this->insert(std::make_pair("C8", std::make_unique<C8_susy>()));  this->insert(std::make_pair("C9", std::make_unique<C9_susy>())); 
        this->insert(std::make_pair("C10", std::make_unique<C10_susy>())); 
    }
    BCoefficientGroup_susy(double Q_match) {
        this->insert(std::make_pair("C1", std::make_unique<C1_susy>(Q_match))); this->insert(std::make_pair("C2", std::make_unique<C2_susy>(Q_match))); this->insert(std::make_pair("C3", std::make_unique<C3_susy>(Q_match)));
        this->insert(std::make_pair("C4", std::make_unique<C4_susy>(Q_match)));  this->insert(std::make_pair("C5", std::make_unique<C5_susy>(Q_match))); this->insert(std::make_pair("C6", std::make_unique<C6_susy>(Q_match))); 
        this->insert(std::make_pair("C7", std::make_unique<C7_susy>(Q_match)));  this->insert(std::make_pair("C8", std::make_unique<C8_susy>(Q_match)));  this->insert(std::make_pair("C9", std::make_unique<C9_susy>(Q_match))); 
        this->insert(std::make_pair("C10", std::make_unique<C10_susy>(Q_match)));
    }
    // BCoefficientGroup_susy(double Q_match, int gen) {
    //     this->insert(std::make_pair("C1", std::make_unique<C1_susy>(Q_match, gen))); this->insert(std::make_pair("C2", std::make_unique<C2_susy>(Q_match, gen))); this->insert(std::make_pair("C3", std::make_unique<C3_susy>(Q_match, gen)));
    //     this->insert(std::make_pair("C4", std::make_unique<C4_susy>(Q_match, gen)));  this->insert(std::make_pair("C5", std::make_unique<C5_susy>(Q_match, gen))); this->insert(std::make_pair("C6", std::make_unique<C6_susy>(Q_match, gen))); 
    //     this->insert(std::make_pair("C7", std::make_unique<C7_susy>(Q_match, gen)));  this->insert(std::make_pair("C8", std::make_unique<C8_susy>(Q_match, gen)));  this->insert(std::make_pair("C9", std::make_unique<C9_susy>(Q_match, gen))); 
    //     this->insert(std::make_pair("C10", std::make_unique<C10_susy>(Q_match, gen)));
    // }

};

class BPrimeCoefficientGroup_susy : public BPrimeCoefficientGroup {
public:
    BPrimeCoefficientGroup_susy() {
        this->insert(std::make_pair("CP1", std::make_unique<CP1_susy>())); this->insert(std::make_pair("CP2", std::make_unique<CP2_susy>())); this->insert(std::make_pair("CP3", std::make_unique<CP3_susy>()));
        this->insert(std::make_pair("CP4", std::make_unique<CP4_susy>()));  this->insert(std::make_pair("CP5", std::make_unique<CP5_susy>())); this->insert(std::make_pair("CP6", std::make_unique<CP6_susy>())); 
        this->insert(std::make_pair("CP7", std::make_unique<CP7_susy>()));  this->insert(std::make_pair("CP8", std::make_unique<CP8_susy>()));  this->insert(std::make_pair("CP9", std::make_unique<CP9_susy>())); 
        this->insert(std::make_pair("CP10", std::make_unique<CP10_susy>())); this->insert(std::make_pair("CPQ1", std::make_unique<CPQ1_susy>())); this->insert(std::make_pair("CPQ2", std::make_unique<CPQ2_susy>())); 
    }
    BPrimeCoefficientGroup_susy(double Q_match) {
        this->insert(std::make_pair("CP1", std::make_unique<CP1_susy>(Q_match))); this->insert(std::make_pair("CP2", std::make_unique<CP2_susy>(Q_match))); this->insert(std::make_pair("CP3", std::make_unique<CP3_susy>(Q_match)));
        this->insert(std::make_pair("CP4", std::make_unique<CP4_susy>(Q_match)));  this->insert(std::make_pair("CP5", std::make_unique<CP5_susy>(Q_match))); this->insert(std::make_pair("CP6", std::make_unique<CP6_susy>(Q_match))); 
        this->insert(std::make_pair("CP7", std::make_unique<CP7_susy>(Q_match)));  this->insert(std::make_pair("CP8", std::make_unique<CP8_susy>(Q_match)));  this->insert(std::make_pair("CP9", std::make_unique<CP9_susy>(Q_match))); 
        this->insert(std::make_pair("CP10", std::make_unique<CP10_susy>(Q_match))); this->insert(std::make_pair("CPQ1", std::make_unique<CPQ1_susy>(Q_match))); this->insert(std::make_pair("CPQ2", std::make_unique<CPQ2_susy>(Q_match)));
    }


};

class BScalarCoefficientGroup_susy : public BScalarCoefficientGroup {
public:
    BScalarCoefficientGroup_susy() : BScalarCoefficientGroup() {
        sus_param = susy_parameters::GetInstance(this->Q_match);
        this->insert(std::make_pair("CQ1", std::make_unique<CQ1_susy>())); this->insert(std::make_pair("CQ2", std::make_unique<CQ2_susy>()));
    }
    BScalarCoefficientGroup_susy(double Q_match) : BScalarCoefficientGroup(Q_match) {
        sus_param = susy_parameters::GetInstance(this->Q_match);
        this->insert(std::make_pair("CQ1", std::make_unique<CQ1_susy>(Q_match))); this->insert(std::make_pair("CQ2", std::make_unique<CQ2_susy>(Q_match)));
    }

    void set_base_1_LO();

private:
    Parameters* sm = Parameters::GetInstance(0);
    Parameters* susy = Parameters::GetInstance(1);
    susy_parameters* sus_param;
};

#endif