#include "Wilsonv2.h"
#include "susy_parameters.h"

class WilsonCoefficient_susy {
protected:
    void set_mod_parameters(Parameters* new_mod) {this->mod = new_mod;};

    WilsonCoefficient_susy (double Q_match) {this->Q_match = Q_match; sus_param = susy_parameters::GetInstance(this->Q_match);
    calculateContribution = [&](auto hFunc, const Array3D_3x7x4& X, const Array3D_3x7x4& X2, int ie, int ae, bool isChargeps) -> double {
            double ratio = std::pow((*sm)("MASS", 24) / (*sus_param).Mch[ie], 2);
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

	susy_parameters* sus_param;;
};
class C1_susy : public C1, public WilsonCoefficient_susy {
public:
    C1_susy(double Q_match) : C1(Q_match), WilsonCoefficient_susy(Q_match) {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation()  {return {0,0};} 
    std::complex<double> NNLO_calculation();

};

class C2_susy : public C2, public WilsonCoefficient_susy {
public:
    C2_susy(double Q_match) : C2(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation()  {return {0,0};} 
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C3_susy : public C3, public WilsonCoefficient_susy {
public:
    C3_susy(double Q_match) : C3(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation() {return {0,0};} 
    std::complex<double> NNLO_calculation();

};

class C4_susy : public C4, public WilsonCoefficient_susy {
public:
    C4_susy(double Q_match) : C4(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation(); 

};

class C5_susy : public C5, public WilsonCoefficient_susy {
public:
    C5_susy(double Q_match) : C5(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation() {return {0,0};} 
    std::complex<double> NNLO_calculation();

};

class C6_susy : public C6, public WilsonCoefficient_susy {
public:
    C6_susy(double Q_match) : C6(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation() {return {0,0};} 
    std::complex<double> NNLO_calculation();

};

class C7_susy : public C7, public WilsonCoefficient_susy {
public:
    C7_susy(double Q_match) : C7(Q_match), WilsonCoefficient_susy(Q_match)  {}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C8_susy : public C8, public WilsonCoefficient_susy {
public:
    C8_susy(double Q_match) : C8(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C9_susy : public C9, public WilsonCoefficient_susy {
public:
    C9_susy(double Q_match) : C9(Q_match), WilsonCoefficient_susy(Q_match)  {}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C10_susy : public C10, public WilsonCoefficient_susy {
public:
    C10_susy(double Q_match) : C10(Q_match), WilsonCoefficient_susy(Q_match)  {}


    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CQ1_susy : public CQ1, public WilsonCoefficient_susy {
public:
    CQ1_susy(double Q_match) : CQ1(Q_match), WilsonCoefficient_susy(Q_match) {}
    CQ1_susy(double Q_match, int gen) : CQ1(Q_match), WilsonCoefficient_susy(Q_match) {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CQ2_susy : public CQ2, public WilsonCoefficient_susy {
public:
    CQ2_susy(double Q_match) : CQ2(Q_match), WilsonCoefficient_susy(Q_match) {}
    CQ2_susy(double Q_match, int gen) : CQ2(Q_match), WilsonCoefficient_susy(Q_match) {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP1_susy : public CP1, public WilsonCoefficient_susy {
public:
    CP1_susy(double Q_match) : CP1(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP1_susy(double Q_match, int gen) : CP1(Q_match), WilsonCoefficient_susy(Q_match) {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP2_susy : public CP2, public WilsonCoefficient_susy {
public:
    CP2_susy(double Q_match) : CP2(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP2_susy(double Q_match, int gen) : CP2(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP3_susy : public CP3, public WilsonCoefficient_susy {
public:
    CP3_susy(double Q_match) : CP3(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP3_susy(double Q_match, int gen) : CP3(Q_match), WilsonCoefficient_susy(Q_match) {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP4_susy : public CP4, public WilsonCoefficient_susy {
public:
    CP4_susy(double Q_match) : CP4(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP4_susy(double Q_match, int gen) : CP4(Q_match), WilsonCoefficient_susy(Q_match) {}


    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP5_susy : public CP5, public WilsonCoefficient_susy {
public:
    CP5_susy(double Q_match) : CP5(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP5_susy(double Q_match, int gen) : CP5(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP6_susy : public CP6, public WilsonCoefficient_susy {
public:
    CP6_susy(double Q_match) : CP6(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP6_susy(double Q_match, int gen) : CP6(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP7_susy : public CP7, public WilsonCoefficient_susy {
public:
    CP7_susy(double Q_match) : CP7(Q_match) , WilsonCoefficient_susy(Q_match){}
    CP7_susy(double Q_match, int gen) : CP7(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP8_susy : public CP8, public WilsonCoefficient_susy {
public:
    CP8_susy(double Q_match) : CP8(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP8_susy(double Q_match, int gen) : CP8(Q_match) , WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP9_susy : public CP9, public WilsonCoefficient_susy {
public:
    CP9_susy(double Q_match) : CP9(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP9_susy(double Q_match, int gen) : CP9(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CP10_susy : public CP10, public WilsonCoefficient_susy {
public:
    CP10_susy(double Q_match) : CP10(Q_match), WilsonCoefficient_susy(Q_match) {}
    CP10_susy(double Q_match, int gen) : CP10(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CPQ1_susy : public CPQ1, public WilsonCoefficient_susy {
public:
    CPQ1_susy(double Q_match) : CPQ1(Q_match), WilsonCoefficient_susy(Q_match) {}
    CPQ1_susy(double Q_match, int gen) : CPQ1(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class CPQ2_susy : public CPQ2, public WilsonCoefficient_susy {
public:
    CPQ2_susy(double Q_match) : CPQ2(Q_match), WilsonCoefficient_susy(Q_match) {}
    CPQ2_susy(double Q_match, int gen) : CPQ2(Q_match), WilsonCoefficient_susy(Q_match) {}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};