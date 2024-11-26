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
        throw std::out_of_range("PDG code not found in " + this->blockname);
    }

    void setValue(int pdgCode, double value, bool force = false) override {
        JSONParser::getInstance(0)->addElement(this->blockname.substr(0, this->blockname.size()-5), pdgCode, value);
        Parameter param (this->blockname.substr(0, this->blockname.size()-5), pdgCode, value, 0);
        if (force) {
            values[pdgCode] = param;
        } else {
            values.emplace(std::make_pair(pdgCode, param));
        }
    }

    void setMode(int pdgCode, ParameterMode mode) override {
        values.at(pdgCode).set_mode(mode);
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
        auto p = Parameter(this->blockname.substr(0, this->blockname.size()-5), pdgCode, value, 0);
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


protected:
    std::array<std::array<Parameter, column>, index> values;
};


// Concrete block for masses
class MassBlock : public MapBlock{
public:
    MassBlock() {this->blockname = "MASSBlock";}
};

// Concrete block for gauge parameters
class GaugeBlock : public MapBlock{
public:
    GaugeBlock() {this->blockname = "GAUGEBlock";}
};

// Concrete block for sminputs parameters
class SMInputBlock : public MapBlock {
public:
    SMInputBlock() {this->blockname = "SMINPUTBlock";}
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

struct WCoef {
    BWilsonCoefficients name;
    complex_t value_LO;
    complex_t value_NLO;
    complex_t value_NNLO;
    double scale;
};

class WilsonBlock : public Block {
public:
    WCoef getValue(BWilsonCoefficients name) const {
        auto it = values.find(name);
        if (it != values.end()) {
            return it->second;
        }
        throw std::out_of_range("PDG code not found in " + this->blockname);
    }

    double getValue(int pdgCode) const {return 0;}
    void setValue(int pdgCode, double value, bool force = false) {}

    void setValue(BWilsonCoefficients name, complex_t value, int order, int type) {
        values[name] = {name, type, order, value};
    }

    void setMode(int pdgCode, ParameterMode mode) {}
protected:
    std::map<BWilsonCoefficients, WCoef> values;

};