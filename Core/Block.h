#pragma once

#include <map>
#include <string>
#include <memory>
#include <stdexcept>
#include <array>

// Abstract base class for all blocks
class Block {
public:
    virtual double getValue(int pdgCode) const = 0;
    virtual void setValue(int pdgCode, double value) = 0;
    virtual ~Block() = default;
    std::string blockname{};
};

class MapBlock : public Block {
public:
    double getValue(int pdgCode) const override {
        auto it = values.find(pdgCode);
        if (it != values.end()) {
            return it->second;
        }
        std::cout << "pdg code is : " << std::endl;
        std::cout << pdgCode << std::endl;
        throw std::out_of_range("PDG code not found in " + this->blockname);
    }

    void setValue(int pdgCode, double value) override {
        values[pdgCode] = value;
    }

protected:
    std::map<int, double> values;

};

template<std::size_t index, std::size_t column>
class ArrayBlock : public Block {
public:
    double getValue(int pdgCode) const override {
        return values[pdgCode/10][pdgCode%10];
    }

    void setValue(int pdgCode, double value) override {
        values[pdgCode/10][pdgCode%10] = value;
    }
    ArrayBlock& operator=(const std::array<std::array<double, column>, index> block) {
        values = block;
        return *this;
    }
protected:
    std::array<std::array<double, column>, index> values;
};


// Concrete block for masses
class MassBlock : public MapBlock{
public:
    MassBlock() {this->blockname = "MassBlock";}
};

// Concrete block for gauge parameters
class GaugeBlock : public MapBlock{
public:
    GaugeBlock() {this->blockname = "GaugeBlock";}
};

// Concrete block for sminputs parameters
class SMInputBlock : public MapBlock {
public:
    SMInputBlock() {this->blockname = "SMInputBlock";}
};

// Concrete block for sminputs parameters
class CKMBlock : public ArrayBlock<3,3> {
public:
    CKMBlock() {this->blockname = "CKMBlock";}
};
/* ------------------------------------------------------------------------------------------------
SUSY BLOCKS*/

class HMIXBlock : public MapBlock {
public:
    HMIXBlock() {this->blockname = "HMIXBlock";}
};

class StopMixBlock : public ArrayBlock<2,2> {
public:
    StopMixBlock() {this->blockname = "StopMixBlock";}
};

class SbotMixBlock : public ArrayBlock<2,2> {
public:
    SbotMixBlock() {this->blockname = "SbotMixBlock";}
};

class StauMixBlock : public ArrayBlock<2,2> {
public:
    StauMixBlock() {this->blockname = "StauMixBlock";}
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
    A0mixBlock() {this->blockname = "A0mixBlock";}
};

class H0mixBlock : public ArrayBlock<4,4> {
public:
    H0mixBlock() {this->blockname = "H0mixBlock";}
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
    AlphaBlock() {this->blockname = "AlphaBlock";}
};

class MSOFTBlock : public MapBlock {
public:
    MSOFTBlock() {this->blockname = "MSOFTBlock";}
};
/* ------------------------------------------------------------------------------------------------
THDM BLOCKS*/

/* ------------------------------------------------------------------------------------------------
Flavor BLOCKS*/

class FlavorBlock {
public:
    double getValue(std::string pdgCode) const {
        auto it = values.find(pdgCode);
        if (it != values.end()) {
            return it->second;
        }
        throw std::out_of_range("PDG code not found in " + this->blockname);
    }

    void setValue(std::string pdgCode, double value) {
        values[pdgCode] = value;
    }

protected:
    std::map<std::string, double> values;
    std::string blockname{};

};

class LifeTimeBlock : public FlavorBlock {
public:
    LifeTimeBlock() {this->blockname = "LifeTimeBlock";}
};

class FConstBlock : public FlavorBlock {
public:
    FConstBlock() {this->blockname = "FConstBlock";}
};

class FLifeBlock : public MapBlock {
public:
    FLifeBlock() {this->blockname = "FLifeBlock";}
};