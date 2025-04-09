#ifndef WILSON_GROUP_H
#define WILSON_GROUP_H

#include <unordered_map>
#include "Include.h"
#include "Utils.h"
#include "Wilson.h"
#include "BWilson.h"
#include "MartyWilson.h"
#include "UseMarty.h"

using BRP = BWilsonRunningParameters;

class CoefficientGroup : public std::map<std::string, std::shared_ptr<WilsonCoefficient>> {
public:
    CoefficientGroup() = default;
    CoefficientGroup(const CoefficientGroup&) = default;
    CoefficientGroup(CoefficientGroup&&) = default;

    CoefficientGroup(std::map<std::string, std::shared_ptr<WilsonCoefficient>>& coeffs) {
        for (auto& coeff : coeffs) {
            this->insert(std::make_pair(coeff.first, std::move(coeff.second)));
        }
        // TODO : retrieve max order of given coeffs and init at this order.
    }

    std::map<ParamId, double> param_cache;

    bool is_double_base() { return this->double_base; }
    virtual void switch_basis() {}
    virtual ~CoefficientGroup() = default;

    void init(QCDOrder order) {
        this->claim_coefficients();
        for (auto& coeff : *this) {
                coeff.second->init(order);
        }
    }

    complex_t getMatching(std::string coeff, std::string order) { return this->at(coeff)->get_CoefficientMatchingValue(order); }
    complex_t getRun(std::string coeff, std::string order) { /* TODO */ }

    void claim_coefficients();

    virtual void init_running_block(QCDOrder order, BWilsonBasis basis = BWilsonBasis::STANDARD) = 0;

    virtual std::shared_ptr<CoefficientGroup> clone() const = 0;

protected:
    bool double_base = false;
};


class BCoefficientGroup : public CoefficientGroup {

public:
    BCoefficientGroup() {
        LOG_INFO("In BCoefficientGroup constructor");
        init_running_parameter_blocks();
        
        if (UseMarty().get()) {
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

    static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>);
    static void base_2_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>);
    static void base_1_NLO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>);
    static void base_2_NLO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>);
    static void base_1_NNLO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>);

    void init_running_block(QCDOrder order, BWilsonBasis basis) override;

    void set_gen(int new_gen) {}

    void switch_basis() override;

    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<BCoefficientGroup>(*this);
    }

protected:
    bool double_base = true;
    QCDOrder current_order = QCDOrder::NONE;
    BWilsonBasis basis = BWilsonBasis::STANDARD;

private:
    void init_running_parameter_blocks();
};


class BPrimeCoefficientGroup : public CoefficientGroup {
public:
    BPrimeCoefficientGroup() {
        if (UseMarty().get()) {
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
    
    static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>);

    void init_running_block(QCDOrder order, BWilsonBasis basis) override;

    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<BPrimeCoefficientGroup>(*this);
    }
};

class BScalarCoefficientGroup : public CoefficientGroup {
public:
    BScalarCoefficientGroup() {
        if (UseMarty().get()) {
            for (auto&& coeff : {"CQ1", "CQ2"}) {
                this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(coeff)));
            }
            return;
        }
        this->insert(std::make_pair("CQ1", std::make_shared<CQ1>())); this->insert(std::make_pair("CQ2", std::make_shared<CQ2>()));
    }

    static void base_1_LO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>);
    static void base_1_NLO_calculation(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>);

    void init_running_block(QCDOrder order, BWilsonBasis basis) override;

    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<BScalarCoefficientGroup>(*this);
    }
};

class BlnuCoefficientGroup : public CoefficientGroup {
public:
    BlnuCoefficientGroup() {
        if (UseMarty().get()) {
            for (auto&& coeff : {"C_Blnu_A", "C_Blnu_P"}) {
                this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(coeff)));
            }
            return;
        }
        this->insert(std::make_pair("C_Blnu_A", std::make_shared<C_Blnu_A>()));
        this->insert(std::make_pair("C_Blnu_P", std::make_shared<C_Blnu_P>()));
    }

    void init_running_block(QCDOrder order, BWilsonBasis basis = BWilsonBasis::STANDARD) {}

    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<BlnuCoefficientGroup>(*this);
    }
};

class BclnuCoefficientGroup : public CoefficientGroup {
    public:
        BclnuCoefficientGroup() {
            if (UseMarty().get()) {
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

        void init_running_block(QCDOrder order, BWilsonBasis basis = BWilsonBasis::STANDARD) {}
    
        std::shared_ptr<CoefficientGroup> clone() const override {
            return std::make_shared<BclnuCoefficientGroup>(*this);
        }
    };

inline std::ostream& operator<<(std::ostream& os, const CoefficientGroup& coeffs) {
    for(auto& coeff : coeffs) {
        os << coeff.second->get_name() << " --------------------------------" << std::endl;
        os << "Matching value at LO: " << coeff.second->get_CoefficientMatchingValue("LO") << std::endl;
        os << "Running value at LO: " << coeff.second->get_CoefficientRunValue("LO") << std::endl;
        os << "Matching value at NLO: " << coeff.second->get_CoefficientMatchingValue("NLO") << std::endl;
        os << "Running value at NLO: " << coeff.second->get_CoefficientRunValue("NLO") << std::endl;
        os << "Matching value at NNLO: " << coeff.second->get_CoefficientMatchingValue("NNLO") << std::endl;
        os << "Running value at NNLO: " << coeff.second->get_CoefficientRunValue("NNLO") << std::endl;
    }
    return os;
}

inline std::ostream& operator<<(std::ostream& os, std::shared_ptr<CoefficientGroup>& coeffs) {
    return os << *coeffs;
}

#endif