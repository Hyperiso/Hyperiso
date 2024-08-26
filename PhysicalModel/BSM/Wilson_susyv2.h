#include "Wilsonv2.h"
#include "susy_parameters.h"

class WilsonCoefficient_susy {
protected:
    void set_mod_parameters(Parameters* new_mod) {this->mod = new_mod;};
    void set_gen(int new_gen) {this->gen = new_gen; susy_params->set_gen(new_gen);}

    Parameters* mod = Parameters::GetInstance(2);
    int gen{3};

    susy_parameters *susy_params = susy_parameters::GetInstance();
};
class C1_susy : public C1, public WilsonCoefficient_susy {
public:
    C1_susy(double Q_match) : C1(Q_match) {}
    C1_susy(double Q_match, int gen) : C1(Q_match) {this->gen = gen;
    susy_params->set_gen(gen);
    }
    C1_susy() : C1() {

    }

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C2_susy : public C2, public WilsonCoefficient_susy {
public:
    C2_susy(double Q_match) : C2(Q_match) {}
    C2_susy(double Q_match, int gen) : C2(Q_match) {this->gen = gen;
    susy_params->set_gen(gen);
    }
    C2_susy() : C2() {

    }

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C3_susy : public C3, public WilsonCoefficient_susy {
public:
    C3_susy(double Q_match) : C3(Q_match) {}
    C3_susy(double Q_match, int gen) : C3(Q_match) {this->gen = gen;
    susy_params->set_gen(gen);
    }
    C3_susy() : C3() {

    }

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C4_susy : public C4, public WilsonCoefficient_susy {
public:
    C4_susy(double Q_match) : C4(Q_match) {}
    C4_susy(double Q_match, int gen) : C4(Q_match) {this->gen = gen;
    susy_params->set_gen(gen);
    }
    C4_susy() : C4() {

    }

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C5_susy : public C5, public WilsonCoefficient_susy {
public:
    C5_susy(double Q_match) : C5(Q_match) {}
    C5_susy(double Q_match, int gen) : C5(Q_match) {this->gen = gen;
    susy_params->set_gen(gen);
    }
    C5_susy() : C5() {

    }

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C6_susy : public C6, public WilsonCoefficient_susy {
public:
    C6_susy(double Q_match) : C6(Q_match) {}
    C6_susy(double Q_match, int gen) : C6(Q_match) {this->gen = gen;
    susy_params->set_gen(gen);
    }
    C6_susy() : C6() {

    }

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C7_susy : public C7, public WilsonCoefficient_susy {
public:
    C7_susy(double Q_match) : C7(Q_match) {}
    C7_susy(double Q_match, int gen) : C7(Q_match) {this->gen = gen;
    susy_params->set_gen(gen);
    }
    C7_susy() : C7() {

    }

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C8_susy : public C8, public WilsonCoefficient_susy {
public:
    C8_susy(double Q_match) : C8(Q_match) {}
    C8_susy(double Q_match, int gen) : C8(Q_match) {this->gen = gen;
    susy_params->set_gen(gen);
    }
    C8_susy() : C8() {

    }

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C9_susy : public C9, public WilsonCoefficient_susy {
public:
    C9_susy(double Q_match) : C9(Q_match) {}
    C9_susy(double Q_match, int gen) : C9(Q_match) {this->gen = gen;
    susy_params->set_gen(gen);
    }
    C9_susy() : C9() {

    }

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};

class C10_susy : public C10, public WilsonCoefficient_susy {
public:
    C10_susy(double Q_match) : C10(Q_match) {}
    C10_susy(double Q_match, int gen) : C10(Q_match) {this->gen = gen;
    susy_params->set_gen(gen);
    }
    C10_susy() : C10() {

    }

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation() {return {0,0};} 

};