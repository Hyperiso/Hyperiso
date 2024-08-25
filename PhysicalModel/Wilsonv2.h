#include <map>
#include <string>
#include <complex>
#include <vector>

#include "Math.h"
#include "Parameters.h"
#include "Logger.h"
#include "Wilson_parameters.h"

class WilsonCoefficient {
protected:
    WilsonCoefficient() {this->set_Q_match(81.);}
    WilsonCoefficient(double Q_match) { this->set_Q_match(Q_match);}

    void set_CoefficientMatchingValue(std::string order, std::complex<double> CoefficientMatchingValue) {
        this->is_now_calculated(order);
        this->CoefficientMatchingValue[order] = CoefficientMatchingValue;
        }

    void is_now_calculated(std::string order) {this->is_calculated[order] = true;}
    complex_t double_to_complex_save(std::string order, double double_temp) {
        complex_t complex_coeff_temp = {double_temp, 0};
        this->set_CoefficientMatchingValue(order, complex_coeff_temp);
        return complex_coeff_temp;
    }

    std::string CoeffName{};
    Wilson_parameters* W_param = Wilson_parameters::GetInstance();
public:

    void set_Q_match(double Q_match) {
        this->Q_match = Q_match;
        this->W_param->SetMuW(Q_match);}
    void set_Q(double Q) {this->Q = Q; this->W_param->SetMu(Q);}
    void set_name(std::string name) {this->CoeffName = name;}
    void set_Wilson_Parameters(Wilson_parameters* W_param) {this->W_param = W_param;}
    void set_WilsonCoeffRun(std::string order, complex_t value) {this->CoefficientRunValue[order] = value;}
    void set_WilsonCoeffMatching(std::string order ,complex_t value) {this->CoefficientMatchingValue[order] = value;}

    std::complex<double> get_CoefficientMatchingValue(std::string order) const {return this->CoefficientMatchingValue.at(order);}
    std::complex<double> get_CoefficientRunValue(std::string order) const {return this->CoefficientRunValue.at(order);}
    double get_Q_match() const {return this->Q_match;}
    double get_Q() const {return this->Q;}
    Wilson_parameters* get_W_params() const {return this->W_param;}
    std::string get_name() const {return this->CoeffName;}

    virtual std::complex<double> LO_calculation() =0;
    virtual std::complex<double> NLO_calculation() = 0;
    virtual std::complex<double> NNLO_calculation() = 0;

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
    double Q{};
    std::map<std::string, std::complex<double>> CoefficientMatchingValue{{"LO", {0.,0.}}, {"NLO", {0.,0.}}, {"NNLO", {0.,0.}}};
    std::map<std::string, std::complex<double>> CoefficientRunValue{{"LO", {0.,0.}}, {"NLO", {0.,0.}}, {"NNLO", {0.,0.}}};

    std::map<std::string, bool> is_calculated{{"LO", false}, {"NLO", false}, {"NNLO", false}};

    
};

class C1 : public WilsonCoefficient {
public:
    C1(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C1");}
    C1() : WilsonCoefficient() {this->set_name("C1");}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation();


    void set_sm_parameters(Parameters* sm) {this->sm = sm;}

    Parameters* sm = Parameters::GetInstance();
};

class C2 : public WilsonCoefficient {
public:
    C2(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C2");}
    C2() : WilsonCoefficient() {this->set_name("C2");}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation() {return {0.,0.};}
    std::complex<double> NNLO_calculation();


    Parameters* sm = Parameters::GetInstance();
};

class C3 : public WilsonCoefficient {
public:
    C3(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C3");}
    C3() : WilsonCoefficient() {this->set_name("C3");}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation() {return {0,0};} 
    std::complex<double> NNLO_calculation() {return {0,0};} 


    Parameters* sm = Parameters::GetInstance();
};

class C4 : public WilsonCoefficient {
public:
    C4(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C4");}
    C4() : WilsonCoefficient() {this->set_name("C4");}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation();


    Parameters* sm = Parameters::GetInstance();
};

class C5 : public WilsonCoefficient {
public:
    C5(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C5");}
    C5() : WilsonCoefficient() {this->set_name("C5");}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation() {return {0,0};} 
    std::complex<double> NNLO_calculation() {return {0,0};} 


    Parameters* sm = Parameters::GetInstance();
};

class C6 : public WilsonCoefficient {
public:
    C6(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C6");}
    C6() : WilsonCoefficient() {this->set_name("C6");}

    std::complex<double> LO_calculation() {return {0,0};} 
    std::complex<double> NLO_calculation() {return {0,0};} 
    std::complex<double> NNLO_calculation() {return {0,0};} 


    Parameters* sm = Parameters::GetInstance();
};

class C7 : public WilsonCoefficient {
public:
    C7(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C7");}
    C7() : WilsonCoefficient() {this->set_name("C7");}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation();

    Parameters* sm = Parameters::GetInstance();
};

class C8 : public WilsonCoefficient {
public:
    C8(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C8");}
    C8() : WilsonCoefficient() {this->set_name("C8");}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation();

    Parameters* sm = Parameters::GetInstance();
};

class C9 : public WilsonCoefficient {
public:
    C9(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C9");}
    C9() : WilsonCoefficient() {this->set_name("C9");}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation();

    Parameters* sm = Parameters::GetInstance();
};

class C10 : public WilsonCoefficient {
public:
    C10(double Q_match) : WilsonCoefficient(Q_match) {this->set_name("C10");}
    C10() : WilsonCoefficient() {this->set_name("C10");}

    std::complex<double> LO_calculation();
    std::complex<double> NLO_calculation();
    std::complex<double> NNLO_calculation();

    Parameters* sm = Parameters::GetInstance();
};

std::ostream& operator<<(std::ostream& os, WilsonCoefficient& coeff) {
    os << "WilsonCoefficient " << coeff.get_name() << "has matching value (" << coeff.get_Q_match() << " GeV) : " << coeff.get_CoefficientMatchingValue("LO") << " at LO" << std::endl;
    os<< ", " << coeff.get_CoefficientMatchingValue("NLO") << " at NLO, " << coeff.get_CoefficientMatchingValue("NNLO") << "at NNLO" << std::endl;
    os << "and " << "has run value (" << coeff.get_Q() << " GeV) : " << coeff.get_CoefficientRunValue("LO") << " at LO" << std::endl;
    os<< ", " << coeff.get_CoefficientRunValue("NLO") << " at NLO, " << coeff.get_CoefficientRunValue("NNLO") << "at NNLO" << std::endl;
    return os;
}

class CoefficientGroup : public std::map<std::string, std::unique_ptr<WilsonCoefficient>> {
public:
    CoefficientGroup() {}
    CoefficientGroup(std::map<std::string, std::unique_ptr<WilsonCoefficient>>& coeffs) {
        for (auto& coeff : coeffs) {
            this->insert(std::make_pair(coeff.first, std::move(coeff.second)));
        }
    }
    virtual ~CoefficientGroup() = default;

    void init_LO() {
        for (auto& coeff : *this) {
            coeff.second->LO_calculation();
        }
    }

    void init_NLO() {
        for (auto& coeff : *this) {
            coeff.second->NLO_calculation();
        }
    }

    void init_NNLO() {
        for (auto& coeff : *this) {
            coeff.second->NNLO_calculation();
        }
    }

    double get_Q_match() {return this->Q_match;}
    double get_Q_run() {return this->Q_run;}

    void set_Q_match(double Q_match) {this->Q_match = Q_match; for (auto& coeff : *this) {coeff.second->set_Q_match(Q_match);}}
    void set_Q_run(double Q_run) {this->Q_run = Q_run; for (auto& coeff : *this) {coeff.second->set_Q(Q_run);}}

    double Q_match{81};
    double Q_run{81};
};



class BCoefficientGroup : public CoefficientGroup {

public:
    BCoefficientGroup() {
        this->insert(std::make_pair("C1", std::make_unique<C1>())); this->insert(std::make_pair("C2", std::make_unique<C2>())); this->insert(std::make_pair("C3", std::make_unique<C3>()));
        this->insert(std::make_pair("C4", std::make_unique<C4>()));  this->insert(std::make_pair("C5", std::make_unique<C5>())); this->insert(std::make_pair("C6", std::make_unique<C6>())); 
        this->insert(std::make_pair("C7", std::make_unique<C7>()));  this->insert(std::make_pair("C8", std::make_unique<C8>()));  this->insert(std::make_pair("C9", std::make_unique<C9>())); 
        this->insert(std::make_pair("C10", std::make_unique<C10>())); 
    }
    BCoefficientGroup(double Q_match) {
        this->insert(std::make_pair("C1", std::make_unique<C1>(Q_match))); this->insert(std::make_pair("C2", std::make_unique<C2>(Q_match))); this->insert(std::make_pair("C3", std::make_unique<C3>(Q_match)));
        this->insert(std::make_pair("C4", std::make_unique<C4>(Q_match)));  this->insert(std::make_pair("C5", std::make_unique<C5>(Q_match))); this->insert(std::make_pair("C6", std::make_unique<C6>(Q_match))); 
        this->insert(std::make_pair("C7", std::make_unique<C7>(Q_match)));  this->insert(std::make_pair("C8", std::make_unique<C8>(Q_match)));  this->insert(std::make_pair("C9", std::make_unique<C9>(Q_match))); 
        this->insert(std::make_pair("C10", std::make_unique<C10>(Q_match)));
    }

    void set_base_1_LO();
    void set_base_2_LO();

    void set_base_1_NLO();
    void set_base_2_NLO();

    void set_base_1_NNLO();
    void set_base_2_NNLO();

    void set_W_params(Wilson_parameters* new_W_param) {this->W_param = new_W_param; for(auto& coeff : *this) {coeff.second->set_Wilson_Parameters(new_W_param);}}
private:
    Wilson_parameters* W_param = Wilson_parameters::GetInstance();

};

class BScalarPrimeCoefficientGroup : public CoefficientGroup {

    BScalarPrimeCoefficientGroup() {

    }
    BScalarPrimeCoefficientGroup(double Q_match) {

    }

    void set_base_1();
    void set_base_2();

};

std::ostream& operator<<(std::ostream& os, BCoefficientGroup& coeffs) {
    for(auto& coeff : coeffs) {
        os << coeff.second->get_name() << " --------------------------------" << std::endl;
        os << "Matching value at LO (" << coeff.second->get_Q_match() << " GeV) : " << coeff.second->get_CoefficientMatchingValue("LO") << std::endl;
        os << "Running value at LO (" << coeff.second->get_Q() << " GeV) : " << coeff.second->get_CoefficientRunValue("LO") << std::endl;
    }
    return os;
}