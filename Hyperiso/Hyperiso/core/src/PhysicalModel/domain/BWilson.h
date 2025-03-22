#ifndef BWILSON_H
#define BWILSON_H

#include "Wilson.h"

class C1 : public WilsonCoefficient {
public:
    C1(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C1");}
    C1() : WilsonCoefficient() {this->set_name("C1");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation();
    complex_t NNLO_calculation();


};

class C2 : public WilsonCoefficient {
public:
    C2(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C2");}
    C2() : WilsonCoefficient() {this->set_name("C2");}

    complex_t LO_calculation();
    complex_t NLO_calculation() {return {0.,0.};}
    complex_t NNLO_calculation();

};

class C3 : public WilsonCoefficient {
public:
    C3(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C3");}
    C3() : WilsonCoefficient() {this->set_name("C3");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation();

};

class C4 : public WilsonCoefficient {
public:
    C4(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C4");}
    C4() : WilsonCoefficient() {this->set_name("C4");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation();
    complex_t NNLO_calculation();

};

class C5 : public WilsonCoefficient {
public:
    C5(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C5");}
    C5() : WilsonCoefficient() {this->set_name("C5");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation();

};

class C6 : public WilsonCoefficient {
public:
    C6(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C6");}
    C6() : WilsonCoefficient() {this->set_name("C6");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation();

};

class C7 : public WilsonCoefficient {
public:
    C7(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C7");}
    C7() : WilsonCoefficient() {this->set_name("C7");}

    complex_t LO_calculation();
    complex_t NLO_calculation();
    complex_t NNLO_calculation();

};

class C8 : public WilsonCoefficient {
public:
    C8(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C8");}
    C8() : WilsonCoefficient() {this->set_name("C8");}

    complex_t LO_calculation();
    complex_t NLO_calculation();
    complex_t NNLO_calculation();

};

class C9 : public WilsonCoefficient {
public:
    C9(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C9");}
    C9() : WilsonCoefficient() {this->set_name("C9");}

    complex_t LO_calculation();
    complex_t NLO_calculation();
    complex_t NNLO_calculation();

};

class C10 : public WilsonCoefficient {
public:
    C10(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C10");}
    C10() : WilsonCoefficient() {this->set_name("C10");}

    complex_t LO_calculation();
    complex_t NLO_calculation();
    complex_t NNLO_calculation();

};


class CQ1 : public WilsonCoefficient {
public:
    CQ1(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CQ1");}
    CQ1() : WilsonCoefficient() {this->set_name("CQ1");}

    complex_t LO_calculation();
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CQ2 : public WilsonCoefficient {
public:
    CQ2(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CQ2");}
    CQ2() : WilsonCoefficient() {this->set_name("CQ2");}

    complex_t LO_calculation();
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CP1 : public WilsonCoefficient {
public:
    CP1(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP1");}
    CP1() : WilsonCoefficient() {this->set_name("CP1");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CP2 : public WilsonCoefficient {
public:
    CP2(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP2");}
    CP2() : WilsonCoefficient() {this->set_name("CP2");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CP3 : public WilsonCoefficient {
public:
    CP3(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP3");}
    CP3() : WilsonCoefficient() {this->set_name("CP3");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CP4 : public WilsonCoefficient {
public:
    CP4(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP4");}
    CP4() : WilsonCoefficient() {this->set_name("CP4");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CP5 : public WilsonCoefficient {
public:
    CP5(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP5");}
    CP5() : WilsonCoefficient() {this->set_name("CP5");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CP6 : public WilsonCoefficient {
public:
    CP6(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP6");}
    CP6() : WilsonCoefficient() {this->set_name("CP6");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CP7 : public WilsonCoefficient {
public:
    CP7(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP7");}
    CP7() : WilsonCoefficient() {this->set_name("CP7");}

    complex_t LO_calculation();
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CP8 : public WilsonCoefficient {
public:
    CP8(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP8");}
    CP8() : WilsonCoefficient() {this->set_name("CP8");}

    complex_t LO_calculation();
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CP9 : public WilsonCoefficient {
public:
    CP9(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP9");}
    CP9() : WilsonCoefficient() {this->set_name("CP9");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CP10 : public WilsonCoefficient {
public:
    CP10(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP10");}
    CP10() : WilsonCoefficient() {this->set_name("CP10");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CPQ1 : public WilsonCoefficient {
public:
    CPQ1(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CPQ1");}
    CPQ1() : WilsonCoefficient() {this->set_name("CPQ1");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

class CPQ2 : public WilsonCoefficient {
public:
    CPQ2(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CPQ2");}
    CPQ2() : WilsonCoefficient() {this->set_name("CPQ2");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    int gen{2};
};

#endif