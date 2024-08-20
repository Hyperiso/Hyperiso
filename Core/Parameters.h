#pragma once

#include "QCDParameters.h"
#include "BlockAccessor.h"
#include "MemoryManager.h"
#include "Interface.h"

#include <memory>
typedef std::complex<double> complex_t; 



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

class FlAVORModelStrategy : public ModelStrategy {
public:
    void initializeParameters(class Parameters& params) override;
};

class Parameters {
public:
    static Parameters* GetInstance(int modelId = 0);

    double operator()(const std::string& block, int pdgCode);

    double alpha_s(double Q);
    double running_mass(double quarkmass, double Q_init, double Q_end, std::string option_massb = "running", std::string option_masst = "pole");

    // Method to allow ModelStrategy to add blocks
    void addBlock(const std::string& name, std::unique_ptr<Block> block) {
        blockAccessor.addBlock(name, std::move(block));
    }
    void setBlockValue(const std::string& name, int pdgCode, double value) {
        blockAccessor.setValue(name, pdgCode, value);
    }

    void setQCDParameters(const QCDParameters&& qcdparams) {QCDRunner = qcdparams;}

    double get_QCD_masse(std::string masstype);
    double getFlavorParam(FlavorParamType type, const std::string& id);
private:
    explicit Parameters(ModelStrategy* modelStrategy);
    static std::map<int, Parameters*> instances;

    QCDParameters QCDRunner;
    BlockAccessor blockAccessor;
    FlavorBlockAccessor flavorblockAccessor;
    
    ModelStrategy* strategy;

    friend class ParametersFactory;
};

class ParametersFactory {
public:
    static Parameters* GetParameters(int modelId);
private:
    static std::map<int, Parameters*> instances;

    static ModelStrategy* createStrategy(int modelId);
};

std::string doubleToString(double value, int precision);