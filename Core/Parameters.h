#pragma once

#include "QCDParameters.h"
#include "BlockAccessor.h"
#include "MemoryManager.h"
#include "Interface.h"
#include "JsonParameters.h"
#include <memory>

typedef std::complex<double> complex_t; 

constexpr int N_PARAM_INSTANCES = 6;

class ModelStrategy {
public:
    virtual void initializeParameters(class Parameters& params) = 0;
    virtual ~ModelStrategy() = default;
};

class SMModelStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

class SUSYModelStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

class THDMModelStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

class FlavorStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

class GeneralModelStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

class WilsonInputStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

class FormFactorStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

class Parameters {
public:
    static Parameters* GetInstance(ParameterType id = ParameterType::SM);
    static ParameterType GetType(const std::string& block, int pdgCode);
    static double Get(ParamId id);

    double operator()(const std::string& block, int pdgCode);

    double alpha_s(double Q);
    double running_mass(double quarkmass, double Q_init, double Q_end, std::string option_massb = "running", std::string option_masst = "pole");

    // Method to allow ModelStrategy to add blocks
    void addBlock(const std::string& name, std::shared_ptr<Block> block) {
        blockAccessor.addBlock(name, block);
    }

    void setBlockValue(const std::string& name, int pdgCode, double value, bool force = false) {
        if (force && (name == "SMINPUTS")) {
            if(pdgCode ==6) {
                QCDRunner.set_mt_pole(value);
                blockAccessor.setValue(name, pdgCode, this->get_QCD_masse("mt_mt"), force);
                return;
            
            } else if (pdgCode ==5) {
                QCDRunner.set_mb_mb(value);
                blockAccessor.setValue(name, pdgCode, value, force);
                return;
            }
        }
        blockAccessor.setValue(name, pdgCode, value, force);
    }

    void setQCDParameters(const QCDParameters&& qcdparams) {QCDRunner = qcdparams;}

    double get_QCD_masse(std::string masstype);

    void changeParameterMode(const ParamId& param_id, ParameterMode new_mode);
    void shiftParameter(const ParamId& param_id, double shift_value);

    static complex_t get_c_CKM_entry(int idx) {
        auto p = Parameters::GetInstance();
        return complex_t((*p)("RECKM", idx), (*p)("IMCKM", idx));
    }


    QCDParameters* QCDaddress() {
        return &this->QCDRunner;
    }
  
    ~Parameters() { LOG_DEBUG("Parameters at ", this); }

private:
    explicit Parameters(ModelStrategy* modelStrategy);
    static std::map<ParameterType, Parameters*> instances;
    std::map<std::pair<std::string, int>, double> originalValuesCache;

    QCDParameters QCDRunner;
    BlockAccessor blockAccessor;

    ModelStrategy* strategy;

    friend class ParametersFactory;
};

class ParametersFactory {
public:
    static Parameters* GetParameters(ParameterType id);
private:
    static std::map<ParameterType, Parameters*> instances;

    static ModelStrategy* createStrategy(ParameterType id);
};

std::string doubleToString(double value, int precision);