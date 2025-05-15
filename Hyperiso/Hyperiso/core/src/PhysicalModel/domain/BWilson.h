// #ifndef BWILSON_H
// #define BWILSON_H

// #include "Wilson.h"

// class C1 : public WilsonCoefficient {
// public:
//     C1() : WilsonCoefficient("C1", GroupMapper::str(WGroup::B) + "_MATCH") {}

//     void LO_calculation() {}
//     void NLO_calculation();
//     void NNLO_calculation();

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<C1>(*this);
//     }

// };

// class C2 : public WilsonCoefficient {
// public:
//     C2() : WilsonCoefficient("C2", GroupMapper::str(WGroup::B) + "_MATCH") {}

//     void LO_calculation();
//     void NLO_calculation() {}
//     void NNLO_calculation();

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<C2>(*this);
//     }

// };

// class C3 : public WilsonCoefficient {
// public:
//     C3() : WilsonCoefficient("C3", GroupMapper::str(WGroup::B) + "_MATCH") {}

//     void LO_calculation() {}
//     void NLO_calculation() {}
//     void NNLO_calculation();

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<C3>(*this);
//     }
// };

// class C4 : public WilsonCoefficient {
// public:
//     C4() : WilsonCoefficient("C4", GroupMapper::str(WGroup::B) + "_MATCH") {}

//     void LO_calculation() {}
//     void NLO_calculation();
//     void NNLO_calculation();

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<C4>(*this);
//     }

// };

// class C5 : public WilsonCoefficient {
// public:
//     C5() : WilsonCoefficient("C5", GroupMapper::str(WGroup::B) + "_MATCH") {}

//     void LO_calculation() {} 
//     void NLO_calculation() {} 
//     void NNLO_calculation();

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<C5>(*this);
//     }

// };

// class C6 : public WilsonCoefficient {
// public:
//     C6() : WilsonCoefficient("C6", GroupMapper::str(WGroup::B) + "_MATCH") {}

//     void LO_calculation() {} 
//     void NLO_calculation() {} 
//     void NNLO_calculation();

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<C6>(*this);
//     }

// };

// class C7 : public WilsonCoefficient {
// public:
//     C7() : WilsonCoefficient("C7", GroupMapper::str(WGroup::B) + "_MATCH") {}

//     void LO_calculation();
//     void NLO_calculation();
//     void NNLO_calculation();

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<C7>(*this);
//     }
// };

// class C8 : public WilsonCoefficient {
// public:
//     C8() : WilsonCoefficient("C8", GroupMapper::str(WGroup::B) + "_MATCH") {}

//     void LO_calculation();
//     void NLO_calculation();
//     void NNLO_calculation();

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<C8>(*this);
//     }
// };

// class C9 : public WilsonCoefficient {
// public:
//     C9() : WilsonCoefficient("C9", GroupMapper::str(WGroup::B) + "_MATCH") {}

//     void LO_calculation();
//     void NLO_calculation();
//     void NNLO_calculation();

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<C9>(*this);
//     }
// };

// class C10 : public WilsonCoefficient {
// public:
//     C10() : WilsonCoefficient("C10", GroupMapper::str(WGroup::B) + "_MATCH") {}

//     void LO_calculation();
//     void NLO_calculation();
//     void NNLO_calculation();

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<C10>(*this);
//     }
// };


// class CQ1 : public WilsonCoefficient {
// public:
//     CQ1() : WilsonCoefficient("CQ1", GroupMapper::str(WGroup::BScalar) + "_MATCH") {}

//     void LO_calculation();
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CQ1>(*this);
//     }

//     int gen{2};
// };

// class CQ2 : public WilsonCoefficient {
// public:
//     CQ2() : WilsonCoefficient("CQ2", GroupMapper::str(WGroup::BScalar) + "_MATCH") {}

//     void LO_calculation();
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CQ2>(*this);
//     }

//     int gen{2};
// };

// class CP1 : public WilsonCoefficient {
// public:
//     CP1() : WilsonCoefficient("CP1", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

//     void LO_calculation() {} 
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CP1>(*this);
//     }

//     int gen{2};
// };

// class CP2 : public WilsonCoefficient {
// public:
//     CP2() : WilsonCoefficient("CP2", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

//     void LO_calculation() {} 
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CP2>(*this);
//     }

//     int gen{2};
// };

// class CP3 : public WilsonCoefficient {
// public:
//     CP3() : WilsonCoefficient("CP3", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

//     void LO_calculation() {} 
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CP3>(*this);
//     }

//     int gen{2};
// };

// class CP4 : public WilsonCoefficient {
// public:
//     CP4() : WilsonCoefficient("CP4", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

//     void LO_calculation() {} 
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CP4>(*this);
//     }

//     int gen{2};
// };

// class CP5 : public WilsonCoefficient {
// public:
//     CP5() : WilsonCoefficient("CP5", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

//     void LO_calculation() {} 
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CP5>(*this);
//     }

//     int gen{2};
// };

// class CP6 : public WilsonCoefficient {
// public:
//     CP6() : WilsonCoefficient("CP6", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

//     void LO_calculation() {} 
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CP6>(*this);
//     }

//     int gen{2};
// };

// class CP7 : public WilsonCoefficient {
// public:
//     CP7() : WilsonCoefficient("CP7", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

//     void LO_calculation();
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CP7>(*this);
//     }

//     int gen{2};
// };

// class CP8 : public WilsonCoefficient {
// public:
//     CP8() : WilsonCoefficient("CP8", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

//     void LO_calculation();
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CP8>(*this);
//     }

//     int gen{2};
// };

// class CP9 : public WilsonCoefficient {
// public:
//     CP9() : WilsonCoefficient("CP9", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

//     void LO_calculation() {} 
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CP9>(*this);
//     }

//     int gen{2};
// };

// class CP10 : public WilsonCoefficient {
// public:
//     CP10() : WilsonCoefficient("CP10", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

//     void LO_calculation() {} 
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CP10>(*this);
//     }

//     int gen{2};
// };

// class CPQ1 : public WilsonCoefficient {
// public:
//     CPQ1() : WilsonCoefficient("CPQ1", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

//     void LO_calculation() {} 
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CPQ1>(*this);
//     }

//     int gen{2};
// };

// class CPQ2 : public WilsonCoefficient {
// public:
//     CPQ2() : WilsonCoefficient("CPQ2", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

//     void LO_calculation() {} 
//     void NLO_calculation() {} 
//     void NNLO_calculation() {} 

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<CPQ2>(*this);
//     }

//     int gen{2};
// };

// #endif