#ifndef BWILSON_H
#define BWILSON_H

#include "Wilson.h"

class C1 : public WilsonCoefficient {
public:
    C1() : WilsonCoefficient("C1", "B_MATCH") {}

    void LO_calculation() {}
    void NLO_calculation();
    void NNLO_calculation();

};

class C2 : public WilsonCoefficient {
public:
    C2() : WilsonCoefficient() {this->set_name("C2");}

    void LO_calculation();
    void NLO_calculation() {}
    void NNLO_calculation();

};

class C3 : public WilsonCoefficient {
public:
    C3() : WilsonCoefficient() {this->set_name("C3");}

    void LO_calculation() {}
    void NLO_calculation() {}
    void NNLO_calculation();

};

class C4 : public WilsonCoefficient {
public:
    C4() : WilsonCoefficient() {this->set_name("C4");}

    void LO_calculation() {}
    void NLO_calculation();
    void NNLO_calculation();

};

class C5 : public WilsonCoefficient {
public:
    C5() : WilsonCoefficient() {this->set_name("C5");}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation();

};

class C6 : public WilsonCoefficient {
public:
    C6() : WilsonCoefficient() {this->set_name("C6");}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation();

};

class C7 : public WilsonCoefficient {
public:
    C7() : WilsonCoefficient() {this->set_name("C7");}

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation();

};

class C8 : public WilsonCoefficient {
public:
    C8() : WilsonCoefficient() {this->set_name("C8");}

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation();

};

class C9 : public WilsonCoefficient {
public:
    C9() : WilsonCoefficient() {this->set_name("C9");}

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation();

};

class C10 : public WilsonCoefficient {
public:
    C10() : WilsonCoefficient() {this->set_name("C10");}

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation();

};


class CQ1 : public WilsonCoefficient {
public:
    CQ1() : WilsonCoefficient() {this->set_name("CQ1");}

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CQ2 : public WilsonCoefficient {
public:
    CQ2() : WilsonCoefficient() {this->set_name("CQ2");}

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP1 : public WilsonCoefficient {
public:
    CP1() : WilsonCoefficient() {this->set_name("CP1");}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP2 : public WilsonCoefficient {
public:
    CP2() : WilsonCoefficient() {this->set_name("CP2");}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP3 : public WilsonCoefficient {
public:
    CP3() : WilsonCoefficient() {this->set_name("CP3");}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP4 : public WilsonCoefficient {
public:
    CP4() : WilsonCoefficient() {this->set_name("CP4");}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP5 : public WilsonCoefficient {
public:
    CP5() : WilsonCoefficient() {this->set_name("CP5");}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP6 : public WilsonCoefficient {
public:
    CP6() : WilsonCoefficient() {this->set_name("CP6");}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP7 : public WilsonCoefficient {
public:
    CP7() : WilsonCoefficient() {this->set_name("CP7");}

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP8 : public WilsonCoefficient {
public:
    CP8() : WilsonCoefficient() {this->set_name("CP8");}

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP9 : public WilsonCoefficient {
public:
    CP9() : WilsonCoefficient() {this->set_name("CP9");}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP10 : public WilsonCoefficient {
public:
    CP10() : WilsonCoefficient() {this->set_name("CP10");}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CPQ1 : public WilsonCoefficient {
public:
    CPQ1() : WilsonCoefficient() {this->set_name("CPQ1");}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CPQ2 : public WilsonCoefficient {
public:
    CPQ2() : WilsonCoefficient() {this->set_name("CPQ2");}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

#endif