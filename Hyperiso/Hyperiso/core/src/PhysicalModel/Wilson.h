#if !defined(HYPERISO_WILSON_H)
#define HYPERISO_WILSON_H
#include <map>
#include <string>
#include <vector>

#include "Math.h"
#include "Utils.h"
#include "Parameters.h"
#include "Logger.h"
#include "Wilson_parameters.h"


class WilsonCoefficient {
protected:
    WilsonCoefficient() {this->set_Q_match(81.);}
    WilsonCoefficient(double Q_match) {this->set_Q_match(Q_match);}


    void is_now_calculated(std::string order) {this->is_calculated[order] = true;}
    complex_t double_to_complex_save(std::string order, double double_temp) {
        complex_t complex_coeff_temp = {double_temp, 0};
        this->set_CoefficientMatchingValue(order, complex_coeff_temp);
        return complex_coeff_temp;
    }

    std::string CoeffName{};
    Wilson_parameters* W_param = Wilson_parameters::GetInstance();
public:
    void set_CoefficientMatchingValue(std::string order, complex_t CoefficientMatchingValue) {
        this->is_now_calculated(order);
        this->CoefficientMatchingValue[order] = CoefficientMatchingValue;
        }

    void set_Q_match(double Q_match) {
        this->Q_match = Q_match;
        this->W_param->SetMuW(Q_match);}
    void set_Q(double Q) {this->Q = Q; this->W_param->SetMu(Q);}
    void set_name(std::string name) {this->CoeffName = name;}
    void set_Wilson_Parameters(Wilson_parameters* W_param) {this->W_param = W_param;}
    void set_WilsonCoeffRun(std::string order, complex_t value) {this->CoefficientRunValue[order] = value;}
    void set_WilsonCoeffMatching(std::string order ,complex_t value) {this->CoefficientMatchingValue[order] = value;}

    complex_t get_CoefficientMatchingValue(std::string order) const {return this->CoefficientMatchingValue.at(order);}
    complex_t get_CoefficientRunValue(std::string order) const {return this->CoefficientRunValue.at(order);}

    complex_t get_CoefficientFullMatchingValue(std::string order) const {
        double fact = QCDHelper::alpha_s(Q_match) / (4 * M_PI);

        if (order == "LO") {
            return this->get_CoefficientMatchingValue("LO");
        }
        else if (order == "NLO") {
            return this->get_CoefficientMatchingValue("LO") + fact *  this->get_CoefficientMatchingValue("NLO");
        }
        else if (order == "NNLO") {
            return this->get_CoefficientMatchingValue("LO") + fact * this->get_CoefficientMatchingValue("NLO") + fact * fact * this->get_CoefficientMatchingValue("NNLO");
        }
        else {
            LOG_ERROR("ValueError", "Order request for wilson getfullmatching is not possible.");
        }
    }

    complex_t get_CoefficientFullRunValue(std::string order) const {
        double fact = QCDHelper::alpha_s(Q) / (4 * M_PI);

        if (order == "LO") {
            return this->get_CoefficientRunValue("LO");
        }
        else if (order == "NLO") {
            return this->get_CoefficientRunValue("LO") + fact * this->get_CoefficientRunValue("NLO");
        }
        else if (order == "NNLO") {
            return this->get_CoefficientRunValue("LO") + fact *  this->get_CoefficientRunValue("NLO") + fact * fact * this->get_CoefficientRunValue("NNLO");
        }
        else {
            LOG_ERROR("ValueError", "Order request for wilson getfullrun is not possible.");
        }
    }

    double get_Q_match() const {return this->Q_match;}
    double get_Q() const {return this->Q;}
    Wilson_parameters* get_W_params() const {return this->W_param;}
    std::string get_name() const {return this->CoeffName;}

    virtual complex_t LO_calculation() =0;
    virtual complex_t NLO_calculation() = 0;
    virtual complex_t NNLO_calculation() = 0;

    bool fill_from_flha();

    bool is_it_calculated(std::string order) {return this->is_calculated[order];}

    bool operator!=(const WilsonCoefficient& other) {
        return this->CoeffName != other.CoeffName;
    }
    WilsonCoefficient& operator+=(const WilsonCoefficient& other) {

        if ((this->get_Q_match() == other.get_Q_match()) && (this->get_Q() == other.get_Q())) {
            this->CoefficientMatchingValue["LO"] += other.get_CoefficientMatchingValue("LO");
            this->CoefficientMatchingValue["NLO"] += other.get_CoefficientMatchingValue("NLO");
            this->CoefficientMatchingValue["NNLO"] += other.get_CoefficientMatchingValue("NNLO");
            this->CoefficientRunValue["LO"] += other.get_CoefficientRunValue("LO");
            this->CoefficientRunValue["NLO"] += other.get_CoefficientRunValue("NLO");
            this->CoefficientRunValue["NNLO"] += other.get_CoefficientRunValue("NNLO");
        }
        else {
            LOG_ERROR("ValueError", "Not the right scale"); 
        }
        return *this;
    }

    virtual ~WilsonCoefficient() = default;
private:
    double Q_match{81};
    double Q{81};
    bool from_lha {false};
    std::map<std::string, complex_t> CoefficientMatchingValue{{"LO", {0.,0.}}, {"NLO", {0.,0.}}, {"NNLO", {0.,0.}}};
    std::map<std::string, complex_t> CoefficientRunValue{{"LO", {0.,0.}}, {"NLO", {0.,0.}}, {"NNLO", {0.,0.}}};

    std::map<std::string, bool> is_calculated{{"LO", false}, {"NLO", false}, {"NNLO", false}};

    
};

class C1 : public WilsonCoefficient {
public:
    C1(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C1");}
    C1() : WilsonCoefficient() {this->set_name("C1");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation();
    complex_t NNLO_calculation();


    void set_sm_parameters(std::shared_ptr<Parameters> sm) {this->sm = sm;}

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C2 : public WilsonCoefficient {
public:
    C2(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C2");}
    C2() : WilsonCoefficient() {this->set_name("C2");}

    complex_t LO_calculation();
    complex_t NLO_calculation() {return {0.,0.};}
    complex_t NNLO_calculation();


    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C3 : public WilsonCoefficient {
public:
    C3(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C3");}
    C3() : WilsonCoefficient() {this->set_name("C3");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation();


    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C4 : public WilsonCoefficient {
public:
    C4(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C4");}
    C4() : WilsonCoefficient() {this->set_name("C4");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation();
    complex_t NNLO_calculation();


    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C5 : public WilsonCoefficient {
public:
    C5(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C5");}
    C5() : WilsonCoefficient() {this->set_name("C5");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation();


    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C6 : public WilsonCoefficient {
public:
    C6(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C6");}
    C6() : WilsonCoefficient() {this->set_name("C6");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation();


    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C7 : public WilsonCoefficient {
public:
    C7(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C7");}
    C7() : WilsonCoefficient() {this->set_name("C7");}

    complex_t LO_calculation();
    complex_t NLO_calculation();
    complex_t NNLO_calculation();

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C8 : public WilsonCoefficient {
public:
    C8(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C8");}
    C8() : WilsonCoefficient() {this->set_name("C8");}

    complex_t LO_calculation();
    complex_t NLO_calculation();
    complex_t NNLO_calculation();

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C9 : public WilsonCoefficient {
public:
    C9(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C9");}
    C9() : WilsonCoefficient() {this->set_name("C9");}

    complex_t LO_calculation();
    complex_t NLO_calculation();
    complex_t NNLO_calculation();

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C10 : public WilsonCoefficient {
public:
    C10(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C10");}
    C10() : WilsonCoefficient() {this->set_name("C10");}

    complex_t LO_calculation();
    complex_t NLO_calculation();
    complex_t NNLO_calculation();

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};


class CQ1 : public WilsonCoefficient {
public:
    CQ1(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CQ1");}
    CQ1() : WilsonCoefficient() {this->set_name("CQ1");}

    complex_t LO_calculation();
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CQ2 : public WilsonCoefficient {
public:
    CQ2(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CQ2");}
    CQ2() : WilsonCoefficient() {this->set_name("CQ2");}

    complex_t LO_calculation();
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CP1 : public WilsonCoefficient {
public:
    CP1(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP1");}
    CP1() : WilsonCoefficient() {this->set_name("CP1");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CP2 : public WilsonCoefficient {
public:
    CP2(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP2");}
    CP2() : WilsonCoefficient() {this->set_name("CP2");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CP3 : public WilsonCoefficient {
public:
    CP3(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP3");}
    CP3() : WilsonCoefficient() {this->set_name("CP3");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CP4 : public WilsonCoefficient {
public:
    CP4(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP4");}
    CP4() : WilsonCoefficient() {this->set_name("CP4");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CP5 : public WilsonCoefficient {
public:
    CP5(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP5");}
    CP5() : WilsonCoefficient() {this->set_name("CP5");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CP6 : public WilsonCoefficient {
public:
    CP6(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP6");}
    CP6() : WilsonCoefficient() {this->set_name("CP6");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CP7 : public WilsonCoefficient {
public:
    CP7(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP7");}
    CP7() : WilsonCoefficient() {this->set_name("CP7");}

    complex_t LO_calculation();
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CP8 : public WilsonCoefficient {
public:
    CP8(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP8");}
    CP8() : WilsonCoefficient() {this->set_name("CP8");}

    complex_t LO_calculation();
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CP9 : public WilsonCoefficient {
public:
    CP9(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP9");}
    CP9() : WilsonCoefficient() {this->set_name("CP9");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CP10 : public WilsonCoefficient {
public:
    CP10(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CP10");}
    CP10() : WilsonCoefficient() {this->set_name("CP10");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CPQ1 : public WilsonCoefficient {
public:
    CPQ1(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CPQ1");}
    CPQ1() : WilsonCoefficient() {this->set_name("CPQ1");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class CPQ2 : public WilsonCoefficient {
public:
    CPQ2(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("CPQ2");}
    CPQ2() : WilsonCoefficient() {this->set_name("CPQ2");}

    complex_t LO_calculation() {return {0,0};} 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
    int gen{2};
};

class C_Blnu_A : public WilsonCoefficient {
public:
    C_Blnu_A(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C_Blnu_A");}
    C_Blnu_A() : WilsonCoefficient() {this->set_name("C_Blnu_A");}

    complex_t LO_calculation() { return double_to_complex_save("LO", 1.); }; 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C_Blnu_P : public WilsonCoefficient {
public:
    C_Blnu_P(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C_Blnu_P");}
    C_Blnu_P() : WilsonCoefficient() {this->set_name("C_Blnu_P");}

    complex_t LO_calculation() {return {0,0};}; 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C_V1 : public WilsonCoefficient {
public:
    C_V1(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C_V1");}
    C_V1() : WilsonCoefficient() {this->set_name("C_V1");}

    complex_t LO_calculation() { return double_to_complex_save("LO", 1.); }; 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C_V2 : public WilsonCoefficient {
public:
    C_V2(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C_V2");}
    C_V2() : WilsonCoefficient() {this->set_name("C_V2");}

    complex_t LO_calculation() {return {0,0};}; 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C_S1 : public WilsonCoefficient {
public:
    C_S1(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C_S1");}
    C_S1() : WilsonCoefficient() {this->set_name("C_S1");}

    complex_t LO_calculation() {return {0,0};}; 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C_S2 : public WilsonCoefficient {
public:
    C_S2(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C_S2");}
    C_S2() : WilsonCoefficient() {this->set_name("C_S2");}

    complex_t LO_calculation() {return {0,0};}; 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

class C_T : public WilsonCoefficient {
public:
    C_T(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C_T");}
    C_T() : WilsonCoefficient() {this->set_name("C_T");}

    complex_t LO_calculation() {return {0,0};}; 
    complex_t NLO_calculation() {return {0,0};} 
    complex_t NNLO_calculation() {return {0,0};} 

    std::shared_ptr<Parameters> sm = Parameters::GetInstance();
};

inline std::ostream& operator<<(std::ostream& os, WilsonCoefficient& coeff) {
    os << "WilsonCoefficient " << coeff.get_name() << "has matching value (" << coeff.get_Q_match() << " GeV) : " << coeff.get_CoefficientMatchingValue("LO") << " at LO" << std::endl;
    os<< ", " << coeff.get_CoefficientMatchingValue("NLO") << " at NLO, " << coeff.get_CoefficientMatchingValue("NNLO") << "at NNLO" << std::endl;
    os << "and " << "has run value (" << coeff.get_Q() << " GeV) : " << coeff.get_CoefficientRunValue("LO") << " at LO" << std::endl;
    os<< ", " << coeff.get_CoefficientRunValue("NLO") << " at NLO, " << coeff.get_CoefficientRunValue("NNLO") << "at NNLO" << std::endl;
    return os;
}



#endif //Wilsonv2