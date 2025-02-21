#ifndef WILSON_GROUP_H
#define WILSON_GROUP_H

#include <map>
#include <string>
#include <complex>
#include <vector>
#include "Wilson.h"
#include "MartyWilson.h"


class CoefficientGroup : public std::map<std::string, std::shared_ptr<WilsonCoefficient>> {
public:

    CoefficientGroup(const CoefficientGroup&) = default;
    CoefficientGroup(CoefficientGroup&&) = default;

    CoefficientGroup() {}
    CoefficientGroup(std::map<std::string, std::shared_ptr<WilsonCoefficient>>& coeffs) {
        for (auto& coeff : coeffs) {
            this->insert(std::make_pair(coeff.first, std::move(coeff.second)));
        }
    }

    std::map<ParamId, double> param_cache;

    bool is_double_base() {return this->double_base;}
    virtual void switch_base() {LOG_ERROR("ValueError", "error");}
    virtual ~CoefficientGroup() = default;


    void init_LO() {
        for (auto& coeff : *this) {
            std::cout << "LO : " << coeff.first << std::endl;
            if (!coeff.second->fill_from_flha())
                coeff.second->LO_calculation();
        }
    }

    void init_NLO() {
        for (auto& coeff : *this) {
            std::cout << "NLO : " << coeff.first << std::endl;
            if (!coeff.second->fill_from_flha())
                coeff.second->NLO_calculation();
        }
    }

    void init_NNLO() {
        for (auto& coeff : *this) {
            std::cout << "NNLO : " << coeff.first << std::endl;
            if (!coeff.second->fill_from_flha())
                coeff.second->NNLO_calculation();
        }
    }

    double get_Q_match() {return this->Q_match;}
    double get_Q_run() {return this->Q_run;}
    std::complex<double> getfullMatching(std::string coeff, std::string order) {return this->at(coeff)->get_CoefficientFullMatchingValue(order);}
    std::complex<double> getfullRun(std::string coeff, std::string order) {return this->at(coeff)->get_CoefficientFullRunValue(order);}

    std::complex<double> getMatching(std::string coeff, std::string order) {return this->at(coeff)->get_CoefficientMatchingValue(order);}
    std::complex<double> getRun(std::string coeff, std::string order) {return this->at(coeff)->get_CoefficientRunValue(order);}

    void setExternalMatchingCoefficient(const std::string& coeff, std::string& order, complex_t value) {
        if (std::shared_ptr<WilsonCoefficient> search = this->find(coeff)->second; search !=this->end()->second) {
            search->set_WilsonCoeffMatching(order, value);
            return;
        }
        LOG_ERROR("KeyError", "matching coefficient", coeff, "Not found in coefficientgroup");
    }

    void setExternalRunningCoefficient(const std::string& coeff, std::string& order, complex_t value) {
        if (std::shared_ptr<WilsonCoefficient> search = this->find(coeff)->second; search !=this->end()->second) {
            search->set_WilsonCoeffRun(order, value);
            return;
        }
        LOG_ERROR("KeyError", "running coefficient", coeff, "Not found in coefficientgroup");
    }

    void set_Q_match(double Q_match) {this->Q_match = Q_match; for (auto& coeff : *this) {coeff.second->set_Q_match(Q_match);}}
    void set_Q_run(double Q_run) {this->Q_run = Q_run; for (auto& coeff : *this) {coeff.second->set_Q(Q_run);}}

    virtual void set_base_1_LO() =0;
    void set_base_2_LO() {}

    virtual void set_base_1_NLO() =0;
    void set_base_2_NLO() {}

    virtual void set_base_1_NNLO() =0;
    void set_base_2_NNLO() {}

    virtual std::shared_ptr<CoefficientGroup> clone() const = 0;
    double Q_match{81};
    double Q_run{81};

    bool double_base = false;

};



class BCoefficientGroup : public CoefficientGroup {

public:
    BCoefficientGroup() {
        if (MemoryManager::GetInstance()->getUseMarty()) {
            for (auto&& coeff : {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"}) {
                this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(coeff)));
            }
            return;
        }
        this->insert(std::make_pair("C1", std::make_shared<C1>())); this->insert(std::make_pair("C2", std::make_shared<C2>())); this->insert(std::make_pair("C3", std::make_shared<C3>()));
        this->insert(std::make_pair("C4", std::make_shared<C4>()));  this->insert(std::make_pair("C5", std::make_shared<C5>())); this->insert(std::make_pair("C6", std::make_shared<C6>())); 
        this->insert(std::make_pair("C7", std::make_shared<C7>()));  this->insert(std::make_pair("C8", std::make_shared<C8>()));  this->insert(std::make_pair("C9", std::make_shared<C9>())); 
        this->insert(std::make_pair("C10", std::make_shared<C10>())); 
    }
    BCoefficientGroup(double Q_match) {
        if (MemoryManager::GetInstance()->getUseMarty()) {
            for (auto&& coeff : {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"}) {
                this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(Q_match, coeff)));
            }
            return;
        }
        this->insert(std::make_pair("C1", std::make_shared<C1>(Q_match))); this->insert(std::make_pair("C2", std::make_shared<C2>(Q_match))); this->insert(std::make_pair("C3", std::make_shared<C3>(Q_match)));
        this->insert(std::make_pair("C4", std::make_shared<C4>(Q_match)));  this->insert(std::make_pair("C5", std::make_shared<C5>(Q_match))); this->insert(std::make_pair("C6", std::make_shared<C6>(Q_match))); 
        this->insert(std::make_pair("C7", std::make_shared<C7>(Q_match)));  this->insert(std::make_pair("C8", std::make_shared<C8>(Q_match)));  this->insert(std::make_pair("C9", std::make_shared<C9>(Q_match))); 
        this->insert(std::make_pair("C10", std::make_shared<C10>(Q_match)));
    }

    void set_base_1_LO();
    void set_base_2_LO();

    void set_base_1_NLO();
    void set_base_2_NLO();

    void set_base_1_NNLO();
    void set_base_2_NNLO();

    void set_W_params(Wilson_parameters* new_W_param) {this->W_param = new_W_param; for(auto& coeff : *this) {coeff.second->set_Wilson_Parameters(new_W_param);}}
    void set_gen(int new_gen) {this->W_param->set_gen(new_gen);}

    void switch_base() override {
        if (base["LO"] == 1) {
            set_base_2_LO();
        }
        else if (base["LO"] == 2) {
            set_base_1_LO();
        }
        if (base["NLO"] == 1) {
            set_base_2_NLO();
        }
        else if (base["NLO"] == 2) {
            set_base_1_NLO();
        }
        if (base["NNLO"] == 1) {
            set_base_2_NNLO();
        }
        else if (base["NNLO"] == 2) {
            set_base_1_NNLO();
        }
    }

    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<BCoefficientGroup>(*this);
    }

protected:
    Wilson_parameters* W_param = Wilson_parameters::GetInstance();
    bool double_base = true;
    std::map<std::string, int> base = {{"LO",0}, {"NLO",0}, {"NNLO",0}};

};


class BPrimeCoefficientGroup : public CoefficientGroup {
public:
    BPrimeCoefficientGroup() {
        if (MemoryManager::GetInstance()->getUseMarty()) {
            for (auto&& coeff : {"CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7", "CP8", "CP9", "CP10", "CPQ1", "CPQ2"}) {
                this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(coeff)));
            }
            return;
        }
        this->insert(std::make_pair("CP1", std::make_shared<CP1>())); this->insert(std::make_pair("CP2", std::make_shared<CP2>())); this->insert(std::make_pair("CP3", std::make_shared<CP3>()));
        this->insert(std::make_pair("CP4", std::make_shared<CP4>()));  this->insert(std::make_pair("CP5", std::make_shared<CP5>())); this->insert(std::make_pair("CP6", std::make_shared<CP6>())); 
        this->insert(std::make_pair("CP7", std::make_shared<CP7>()));  this->insert(std::make_pair("CP8", std::make_shared<CP8>()));  this->insert(std::make_pair("CP9", std::make_shared<CP9>())); 
        this->insert(std::make_pair("CP10", std::make_shared<CP10>())); this->insert(std::make_pair("CPQ1", std::make_shared<CPQ1>())); this->insert(std::make_pair("CPQ2", std::make_shared<CPQ2>())); 
    }
    BPrimeCoefficientGroup(double Q_match) {
        if (MemoryManager::GetInstance()->getUseMarty()) {
            for (auto&& coeff : {"CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7", "CP8", "CP9", "CP10", "CPQ1", "CPQ2"}) {
                this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(Q_match, coeff)));
            }
            return;
        }
        this->insert(std::make_pair("CP1", std::make_shared<CP1>(Q_match))); this->insert(std::make_pair("CP2", std::make_shared<CP2>(Q_match))); this->insert(std::make_pair("CP3", std::make_shared<CP3>(Q_match)));
        this->insert(std::make_pair("CP4", std::make_shared<CP4>(Q_match)));  this->insert(std::make_pair("CP5", std::make_shared<CP5>(Q_match))); this->insert(std::make_pair("CP6", std::make_shared<CP6>(Q_match))); 
        this->insert(std::make_pair("CP7", std::make_shared<CP7>(Q_match)));  this->insert(std::make_pair("CP8", std::make_shared<CP8>(Q_match)));  this->insert(std::make_pair("CP9", std::make_shared<CP9>(Q_match))); 
        this->insert(std::make_pair("CP10", std::make_shared<CP10>(Q_match))); this->insert(std::make_pair("CPQ1", std::make_shared<CPQ1>(Q_match))); this->insert(std::make_pair("CPQ2", std::make_shared<CPQ2>(Q_match)));
    }

    

    void set_base_1_LO();
    void set_base_1_NLO() {}
    void set_base_1_NNLO() {}

    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<BPrimeCoefficientGroup>(*this);
    }

protected:
    Wilson_parameters* W_param = Wilson_parameters::GetInstance();
};

class BScalarCoefficientGroup : public CoefficientGroup {
public:
    BScalarCoefficientGroup() {
        if (MemoryManager::GetInstance()->getUseMarty()) {
            for (auto&& coeff : {"CQ1", "CQ2"}) {
                this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(coeff)));
            }
            return;
        }
        this->insert(std::make_pair("CQ1", std::make_shared<CQ1>())); this->insert(std::make_pair("CQ2", std::make_shared<CQ2>()));
    }
    BScalarCoefficientGroup(double Q_match) {
        if (MemoryManager::GetInstance()->getUseMarty()) {
            for (auto&& coeff : {"CQ1", "CQ2"}) {
                this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(Q_match, coeff)));
            }
            return;
        }
        this->insert(std::make_pair("CQ1", std::make_shared<CQ1>(Q_match))); this->insert(std::make_pair("CQ2", std::make_shared<CQ2>(Q_match)));
    }

    void set_base_1_LO();
    void set_base_1_NLO();
    void set_base_1_NNLO() {}

    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<BScalarCoefficientGroup>(*this);
    }

protected:
    Wilson_parameters* W_param = Wilson_parameters::GetInstance();
};

class BlnuCoefficientGroup : public CoefficientGroup {
public:
    BlnuCoefficientGroup() {
        if (MemoryManager::GetInstance()->getUseMarty()) {
            for (auto&& coeff : {"C_Blnu_A", "C_Blnu_P"}) {
                this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(coeff)));
            }
            return;
        }
        this->insert(std::make_pair("C_Blnu_A", std::make_shared<C_Blnu_A>()));
        this->insert(std::make_pair("C_Blnu_P", std::make_shared<C_Blnu_P>()));
    }

    BlnuCoefficientGroup(double Q_match) {
        if (MemoryManager::GetInstance()->getUseMarty()) {
            for (auto&& coeff : {"C_Blnu_A", "C_Blnu_P"}) {
                this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(Q_match, coeff)));
            }
            return;
        }
        this->insert(std::make_pair("C_Blnu_A", std::make_shared<C_Blnu_A>(Q_match)));
        this->insert(std::make_pair("C_Blnu_P", std::make_shared<C_Blnu_P>(Q_match)));
    }

    void set_base_1_LO() {}
    void set_base_1_NLO() {}
    void set_base_1_NNLO() {}

    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<BlnuCoefficientGroup>(*this);
    }

protected:
    Wilson_parameters* W_param = Wilson_parameters::GetInstance();
};

class BclnuCoefficientGroup : public CoefficientGroup {
    public:
        BclnuCoefficientGroup() {
            if (MemoryManager::GetInstance()->getUseMarty()) {
                for (auto&& coeff : {"C_V1", "C_V2", "C_S1", "C_S2", "C_T"}) {
                    this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(coeff)));
                }
                return;
            }
            this->insert(std::make_pair("C_V1", std::make_shared<C_V1>()));
            this->insert(std::make_pair("C_V2", std::make_shared<C_V2>()));
            this->insert(std::make_pair("C_S1", std::make_shared<C_S1>()));
            this->insert(std::make_pair("C_S2", std::make_shared<C_S2>()));
            this->insert(std::make_pair("C_T", std::make_shared<C_T>()));
        }
    
        BclnuCoefficientGroup(double Q_match) {
            if (MemoryManager::GetInstance()->getUseMarty()) {
                for (auto&& coeff : {"C_V1", "C_V2", "C_S1", "C_S2", "C_T"}) {
                    this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(Q_match, coeff)));
                }
                return;
            }
            this->insert(std::make_pair("C_V1", std::make_shared<C_V1>(Q_match)));
            this->insert(std::make_pair("C_V2", std::make_shared<C_V2>(Q_match)));
            this->insert(std::make_pair("C_S1", std::make_shared<C_S1>(Q_match)));
            this->insert(std::make_pair("C_S2", std::make_shared<C_S2>(Q_match)));
            this->insert(std::make_pair("C_T", std::make_shared<C_T>(Q_match)));
        }
    
        void set_base_1_LO() {}
        void set_base_1_NLO() {}
        void set_base_1_NNLO() {}
    
        std::shared_ptr<CoefficientGroup> clone() const override {
            return std::make_shared<BclnuCoefficientGroup>(*this);
        }
    
    protected:
        Wilson_parameters* W_param = Wilson_parameters::GetInstance();
    };

inline std::ostream& operator<<(std::ostream& os, BCoefficientGroup& coeffs) {
    for(auto& coeff : coeffs) {
        os << coeff.second->get_name() << " --------------------------------" << std::endl;
        os << "Matching value at LO (" << coeff.second->get_Q_match() << " GeV) : " << coeff.second->get_CoefficientMatchingValue("LO") << std::endl;
        os << "Running value at LO (" << coeff.second->get_Q() << " GeV) : " << coeff.second->get_CoefficientRunValue("LO") << std::endl;
        os << "Matching value at NLO (" << coeff.second->get_Q_match() << " GeV) : " << coeff.second->get_CoefficientMatchingValue("NLO") << std::endl;
        os << "Running value at NLO (" << coeff.second->get_Q() << " GeV) : " << coeff.second->get_CoefficientRunValue("NLO") << std::endl;
        os << "Matching value at NNLO (" << coeff.second->get_Q_match() << " GeV) : " << coeff.second->get_CoefficientMatchingValue("NNLO") << std::endl;
        os << "Running value at NNLO (" << coeff.second->get_Q() << " GeV) : " << coeff.second->get_CoefficientRunValue("NNLO") << std::endl;
    }
    return os;
}

inline std::ostream& operator<<(std::ostream& os, std::shared_ptr<CoefficientGroup>& coeffs) {
    for(auto& coeff : *coeffs) {
        os << coeff.second->get_name() << " --------------------------------" << std::endl;
        os << "Matching value at LO (" << coeff.second->get_Q_match() << " GeV) : " << coeff.second->get_CoefficientMatchingValue("LO") << std::endl;
        os << "Running value at LO (" << coeff.second->get_Q() << " GeV) : " << coeff.second->get_CoefficientRunValue("LO") << std::endl;
        os << "Matching value at NLO (" << coeff.second->get_Q_match() << " GeV) : " << coeff.second->get_CoefficientMatchingValue("NLO") << std::endl;
        os << "Running value at NLO (" << coeff.second->get_Q() << " GeV) : " << coeff.second->get_CoefficientRunValue("NLO") << std::endl;
        os << "Matching value at NNLO (" << coeff.second->get_Q_match() << " GeV) : " << coeff.second->get_CoefficientMatchingValue("NNLO") << std::endl;
        os << "Running value at NNLO (" << coeff.second->get_Q() << " GeV) : " << coeff.second->get_CoefficientRunValue("NNLO") << std::endl;
    }
    return os;
}

#endif