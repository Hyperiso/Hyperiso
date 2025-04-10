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
    C2() : WilsonCoefficient("C2", "B_MATCH") {}

    void LO_calculation();
    void NLO_calculation() {}
    void NNLO_calculation();

};

class C3 : public WilsonCoefficient {
public:
    C3() : WilsonCoefficient("C3", "B_MATCH") {}

    void LO_calculation() {}
    void NLO_calculation() {}
    void NNLO_calculation();

};

class C4 : public WilsonCoefficient {
public:
    C4() : WilsonCoefficient("C4", "B_MATCH") {}

    void LO_calculation() {}
    void NLO_calculation();
    void NNLO_calculation();

};

class C5 : public WilsonCoefficient {
public:
    C5() : WilsonCoefficient("C5", "B_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation();

};

class C6 : public WilsonCoefficient {
public:
    C6() : WilsonCoefficient("C6", "B_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation();

};

class C7 : public WilsonCoefficient {
public:
    C7() : WilsonCoefficient("C7", "B_MATCH") {}

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation();

};

class C8 : public WilsonCoefficient {
public:
    C8() : WilsonCoefficient("C8", "B_MATCH") {}

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation();

};

class C9 : public WilsonCoefficient {
public:
    C9() : WilsonCoefficient("C9", "B_MATCH") {}

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation();

};

class C10 : public WilsonCoefficient {
public:
    C10() : WilsonCoefficient("C10", "B_MATCH") {}

    void LO_calculation();
    void NLO_calculation();
    void NNLO_calculation();

};


class CQ1 : public WilsonCoefficient {
public:
    CQ1() : WilsonCoefficient("CQ1", "B_SCALAR_MATCH") {}

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CQ2 : public WilsonCoefficient {
public:
    CQ2() : WilsonCoefficient("CQ2", "B_SCALAR_MATCH") {}

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP1 : public WilsonCoefficient {
public:
    CP1() : WilsonCoefficient("CP1", "B_PRIME_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP2 : public WilsonCoefficient {
public:
    CP2() : WilsonCoefficient("CP2", "B_PRIME_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP3 : public WilsonCoefficient {
public:
    CP3() : WilsonCoefficient("CP3", "B_PRIME_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP4 : public WilsonCoefficient {
public:
    CP4() : WilsonCoefficient("CP4", "B_PRIME_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP5 : public WilsonCoefficient {
public:
    CP5() : WilsonCoefficient("CP5", "B_PRIME_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP6 : public WilsonCoefficient {
public:
    CP6() : WilsonCoefficient("CP6", "B_PRIME_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP7 : public WilsonCoefficient {
public:
    CP7() : WilsonCoefficient("CP7", "B_PRIME_MATCH") {}

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP8 : public WilsonCoefficient {
public:
    CP8() : WilsonCoefficient("CP8", "B_PRIME_MATCH") {}

    void LO_calculation();
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP9 : public WilsonCoefficient {
public:
    CP9() : WilsonCoefficient("CP9", "B_PRIME_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CP10 : public WilsonCoefficient {
public:
    CP10() : WilsonCoefficient("CP10", "B_PRIME_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CPQ1 : public WilsonCoefficient {
public:
    CPQ1() : WilsonCoefficient("CPQ1", "B_PRIME_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

class CPQ2 : public WilsonCoefficient {
public:
    CPQ2() : WilsonCoefficient("CPQ2", "B_PRIME_MATCH") {}

    void LO_calculation() {} 
    void NLO_calculation() {} 
    void NNLO_calculation() {} 

    int gen{2};
};

#endif