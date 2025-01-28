#pragma once

#include <map>
#include <string>
#include <memory>
#include <stdexcept>
#include <array>
#include "JsonParameters.h"
#include "Parameter.h"

// Abstract base class for all blocks
class Block {
public:
    virtual double getValue(int pdgCode) const = 0;
    virtual void setValue(int pdgCode, double value, bool force = false) = 0;
    virtual void setMode(int pdgCode, ParameterMode mode) = 0;
    virtual std::map<int, double> getAllValues() = 0;
    virtual ~Block() = default;
    std::string blockname{};
};

class MapBlock : public Block {
public:
    double getValue(int pdgCode) const override {
        auto it = values.find(pdgCode);
        if (it != values.end()) {
            return it->second.get_val();
        }
        throw std::invalid_argument("PDG code not found in " + this->blockname);
    }

    void setValue(int pdgCode, double value, bool force = false) override {
        JSONParser::getInstance(0)->addElement(this->blockname.substr(0, this->blockname.size()-5), pdgCode, value);
        Parameter param (ParamId {ParameterType::CUSTOM, this->blockname.substr(0, this->blockname.size()-5), pdgCode}, value, 0);
        if (force) {
            values[pdgCode] = param;
        } else {
            // values.emplace(std::make_pair(pdgCode, param));
            values[pdgCode] = param;
        }
    }

    void setMode(int pdgCode, ParameterMode mode) override {
        values.at(pdgCode).set_mode(mode);
    }

    std::map<int, double> getAllValues() override {
        std::map<int, double> map_values;
        for (auto& value : values) {
            map_values[value.first] = value.second.get_val();
        }
        return map_values;
    }

protected:
    std::map<int, Parameter> values;

};

template<std::size_t index, std::size_t column>
class ArrayBlock : public Block {
public:
    double getValue(int pdgCode) const override {
        return values[pdgCode / 10][pdgCode % 10].get_val();
    }

    void setValue(int pdgCode, double value, bool force = false) override {
        JSONParser::getInstance(0)->addElement(this->blockname.substr(0, this->blockname.size()-5), pdgCode, value);
        auto p = Parameter(ParamId {ParameterType::CUSTOM, this->blockname.substr(0, this->blockname.size()-5), pdgCode}, value, 0);
        values[pdgCode / 10][pdgCode % 10] = p;
    }

    void setValues(const std::array<std::array<double, column>, index>& values) {
        for (size_t i = 0; i < index; i++) {
            for (size_t j = 0; j < column; j++) {
                setValue(i * 10 + j, values[i][j]);       
            }
        }
    }

    void setMode(int pdgCode, ParameterMode mode) override {
        values.at(pdgCode/10).at(pdgCode%10).set_mode(mode);
    }

    ArrayBlock& operator=(const std::array<std::array<Parameter, column>, index> block) {
        values = block;
        return *this;
    }

    std::map<int, double> getAllValues() override {
        std::map<int, double> map_values;
        int i,j=0;
        for (auto& value : values) {
            for (auto& valu : value) {
                std::cout << i/10 +j%10 << " " << valu.get_val() << std::endl;
                map_values[i/10 + j++%10] = valu.get_val();
            }
            ++i;
        }
        return map_values;
    }

protected:
    std::array<std::array<Parameter, column>, index> values;
};


// Concrete block for masses
class MassBlock : public MapBlock {
public:
    MassBlock() {this->blockname = "MASSBlock";}

    double getValue(int pdgCode) const override {
        if (pdgCode == 5 || pdgCode == 6) {
            LOG_WARN("Accessing heavy quark masses through Parameters is deprecated. Use QCDHelper instead.");
        }
        return MapBlock::getValue(pdgCode);
    }
};

// Concrete block for gauge parameters
class GaugeBlock : public MapBlock{
public:
    GaugeBlock() {this->blockname = "GAUGEBlock";}
};

// Concrete block for sminputs parameters
class SMInputBlock : public MapBlock {
public:
    SMInputBlock() {this->blockname = "SMINPUTSBlock";}
};

// Concrete block for sminputs parameters
class RECKMBlock : public ArrayBlock<3,3> {
public:
    RECKMBlock() {this->blockname = "RECKMBlock";}
};

class IMCKMBlock : public ArrayBlock<3,3> {
public:
    IMCKMBlock() {this->blockname = "IMCKMBlock";}
};
/* ------------------------------------------------------------------------------------------------
SUSY BLOCKS*/

class HMIXBlock : public MapBlock {
public:
    HMIXBlock() {this->blockname = "HMIXBlock";}
};

class StopMixBlock : public ArrayBlock<2,2> {
public:
    StopMixBlock() {this->blockname = "STOPMIXBlock";}
};

class SbotMixBlock : public ArrayBlock<2,2> {
public:
    SbotMixBlock() {this->blockname = "SBOTMIXBlock";}
};

class StauMixBlock : public ArrayBlock<2,2> {
public:
    StauMixBlock() {this->blockname = "STAUMIXBlock";}
};

class UMIXBlock : public ArrayBlock<2,2> {
public:
    UMIXBlock() {this->blockname = "UMIXBlock";}
};

class VMIXBlock : public ArrayBlock<2,2> {
public:
    VMIXBlock() {this->blockname = "VMIXBlock";}
};

class NMIXBlock : public ArrayBlock<4,4> {
public:
    NMIXBlock() {this->blockname = "NMIXBlock";}
};

class A0mixBlock : public ArrayBlock<4,4> {
public:
    A0mixBlock() {this->blockname = "A0MIXBlock";}
};

class H0mixBlock : public ArrayBlock<4,4> {
public:
    H0mixBlock() {this->blockname = "H0MIXBlock";}
};

class YUBlock : public ArrayBlock<3,3> {
public:
    YUBlock() {this->blockname = "YUBlock";}
};

class YDBlock : public ArrayBlock<3,3> {
public:
    YDBlock() {this->blockname = "YDBlock";}
};

class YEBlock : public ArrayBlock<3,3> {
public:
    YEBlock() {this->blockname = "YEBlock";}
};

class AUBlock : public ArrayBlock<3,3> {
public:
    AUBlock() {this->blockname = "AUBlock";}
};

class ADBlock : public ArrayBlock<3,3> {
public:
    ADBlock() {this->blockname = "ADBlock";}
};

class AEBlock : public ArrayBlock<3,3> {
public:
    AEBlock() {this->blockname = "AEBlock";}
};

class AlphaBlock : public ArrayBlock<1,1> {
public:
    AlphaBlock() {this->blockname = "ALPHABlock";}
};

class MSOFTBlock : public MapBlock {
public:
    MSOFTBlock() {this->blockname = "MSOFTBlock";}
};
/* ------------------------------------------------------------------------------------------------
THDM BLOCKS*/

/* ------------------------------------------------------------------------------------------------
Flavor BLOCKS*/

class FConstBlock : public MapBlock {
public:
    FConstBlock() {this->blockname = "FConstBlock";}
};

class FConstRatioBlock : public MapBlock {
public:
    FConstRatioBlock() {this->blockname = "FConstRatioBlock";}
};

class FBag : public MapBlock {
public:
    FBag() {this->blockname = "FBagBlock";}
};

class FLifeBlock : public MapBlock {
public:
    FLifeBlock() {this->blockname = "FLifeBlock";}
};

/* ------------------------------------------------------------------------------------------------
Wilson Input BLOCKS*/

class WilsonBlock : public Block {
    // pdgCode is NNO with NN = integer ID of the coefficient as in BWilsonCoefficients enum, O is QCD order
    // pdgCode -1 is reserved to access the scale of the coefficients
    // pdgCode -2 is reserved to access the type of the coefficients
public:
    double getValue(int pdgCode) const override {
        if (pdgCode == -1) {
            return scale;
        } else if (pdgCode == -2) {
            return type;
        }

        int order = pdgCode % 10; 
        WCoef id = static_cast<WCoef>((pdgCode - order) / 10); 

        return values.at(id)[order];
    }

    void setValue(int pdgCode, double value, bool force = false) {
        if (pdgCode == -1) {
            scale = value;
        } else if (pdgCode == -2) {
            type = (int)value;
        }

        int order = pdgCode % 10; 
        WCoef id = static_cast<WCoef>((pdgCode - order) / 10); 
        if (!values.contains(id)) {
            values.emplace(std::make_pair(id, std::array<double, 3>()));
        }
        values.at(id)[order] = value;
    }

    void setMode(int pdgCode, ParameterMode mode) {}

    std::map<int, double> getAllValues() override {
        return {};
    }
protected:
    // Index is QCD order
    std::map<WCoef, std::array<double, 3>> values; 
    double scale;
    int type;

};

/* ------------------------------------------------------------------------------------------------
Form Factors BLOCKS*/

class BKsBlock : public MapBlock {
public:
    BKsBlock() {this->blockname = "BKsBlock";}
};

class BllBlock : public MapBlock {
public:
    BllBlock() {this->blockname = "BllBlock";}
};

class BXsBlock : public MapBlock {
public:
    BXsBlock() {this->blockname = "BXsBlock";}
};