/**
 * @file Blocks.h
 * @brief Defines various block classes for parameter management.
 * 
 * This file declares a hierarchy of block classes used to store and manage parameter values in the differents 
 * Parameters instances.
 */
#if !defined(BLOCK_H)
#define BLOCK_H

#include <map>
#include <string>
#include <memory>
#include <stdexcept>
#include <array>
#include "JsonParameters.h"
#include "Parameter.h"

/**
 * @class Block
 * @brief Abstract base class for parameter blocks.
 */
class Block {
public:
    /**
     * @brief Retrieves the value associated with a given PDG code (int).
     * @param pdgCode The PDG code of the parameter.
     * @return The parameter value.
     */
    virtual double getValue(int pdgCode) const = 0;

    /**
     * @brief Sets the value of a parameter.
     * @param pdgCode The PDG code of the parameter.
     * @param value The new value to set.
     * @param force If true, forces the update.
     */
    virtual void setValue(int pdgCode, double value, bool force = false) = 0;

    /**
     * @brief Sets the mode of a parameter.
     * @param pdgCode The PDG code of the parameter.
     * @param mode The mode to set.
     */
    virtual void setMode(int pdgCode, ParameterMode mode) = 0;

    /**
     * @brief Retrieves all parameter values.
     * @return A map of PDG codes to parameter values.
     */
    virtual std::map<int, double> getAllValues() = 0;

    /**
     * @brief Virtual destructor.
     */
    virtual ~Block() = default;

    /// Name of the block
    std::string blockname{};
};

/**
 * @class MapBlock
 * @brief A class for  storing blocks that have only 1 id.
 */
class MapBlock : public Block {
public:
    
    double getValue(int pdgCode) const override;
    void setValue(int pdgCode, double value, bool force = false) override;
    void setMode(int pdgCode, ParameterMode mode) override;
    std::map<int, double> getAllValues() override;

protected:
    /// Map of PDG codes to parameters
    std::map<int, Parameter> values;

};

/**
 * @class ArrayBlock
 * @brief A block storing parameters in a 2D array (for matrix).
 * @tparam index Number of rows.
 * @tparam column Number of columns.
 */
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

    /**
     * @brief Sets values from a 2D array.
     * @param values The array of values to set.
     */
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
        size_t i{0}, j{0};
        for (auto& value : values) {
            for (auto& valu : value) {
                map_values[i * 10 + j++%3] = valu.get_val();
            }
            ++i;
        }
        return map_values;
    }

protected:
    /// 2D array of parameter values
    std::array<std::array<Parameter, column>, index> values;
};

// Specific Block Implementations
/** @class MassBlock @brief Block for mass parameters. */
class MassBlock : public MapBlock {
public:
    MassBlock() {this->blockname = "MASSBlock";}

    double getValue(int pdgCode) const override;
};

/** @class GaugeBlock @brief Block for gauge parameters. */
class GaugeBlock : public MapBlock{
public:
    GaugeBlock() {this->blockname = "GAUGEBlock";}
};

/** @class SMInputBlock @brief Block for Standard Model input parameters. */
class SMInputBlock : public MapBlock {
public:
    SMInputBlock() {this->blockname = "SMINPUTSBlock";}
};

/** @class RECKMBlock @brief Block for real CKM matrix parameters. */
class RECKMBlock : public ArrayBlock<3,3> {
public:
    RECKMBlock() {this->blockname = "RECKMBlock";}
};

/** @class IMCKMBlock @brief Block for imaginary CKM matrix parameters. */
class IMCKMBlock : public ArrayBlock<3,3> {
public:
    IMCKMBlock() {this->blockname = "IMCKMBlock";}
};
/* ------------------------------------------------------------------------------------------------
SUSY BLOCKS*/

/** @class HMIXBlock @brief Block for HMIX parameters in SUSY models */
class HMIXBlock : public MapBlock {
public:
    HMIXBlock() {this->blockname = "HMIXBlock";}
};

/** @class StopMixBlock @brief Block for stop mixing matrix in SUSY models */
class StopMixBlock : public ArrayBlock<2,2> {
public:
    StopMixBlock() {this->blockname = "STOPMIXBlock";}
};

/** @class SbotMixBlock @brief Block for sbot mixing matrix in SUSY models */
class SbotMixBlock : public ArrayBlock<2,2> {
public:
    SbotMixBlock() {this->blockname = "SBOTMIXBlock";}
};

/** @class StauMixBlock @brief Block for stau mixing matrix in SUSY models */
class StauMixBlock : public ArrayBlock<2,2> {
public:
    StauMixBlock() {this->blockname = "STAUMIXBlock";}
};

/** @class UMIXBlock @brief Block for chargino mixing matrix in SUSY models */
class UMIXBlock : public ArrayBlock<2,2> {
public:
    UMIXBlock() {this->blockname = "UMIXBlock";}
};

/** @class VMIXBlock @brief Block for chargino mixing matrix in SUSY models */
class VMIXBlock : public ArrayBlock<2,2> {
public:
    VMIXBlock() {this->blockname = "VMIXBlock";}
};

/** @class NMIXBlock @brief Block for neutralino mixing matrix in SUSY models */
class NMIXBlock : public ArrayBlock<4,4> {
public:
    NMIXBlock() {this->blockname = "NMIXBlock";}
};

/** @class A0mixBlock @brief Block for A0 mixing matrix in SUSY models */
class A0mixBlock : public ArrayBlock<4,4> {
public:
    A0mixBlock() {this->blockname = "A0MIXBlock";}
};

/** @class H0mixBlock @brief Block for H0 mixing matrix in SUSY models */
class H0mixBlock : public ArrayBlock<4,4> {
public:
    H0mixBlock() {this->blockname = "H0MIXBlock";}
};

/** @class YUBlock @brief Block for U Yukawa matrix in SUSY models */
class YUBlock : public ArrayBlock<3,3> {
public:
    YUBlock() {this->blockname = "YUBlock";}
};

/** @class YDBlock @brief Block for D Yukawa matrix in SUSY models */
class YDBlock : public ArrayBlock<3,3> {
public:
    YDBlock() {this->blockname = "YDBlock";}
};

/** @class YEBlock @brief Block for E Yukawa matrix in SUSY models */
class YEBlock : public ArrayBlock<3,3> {
public:
    YEBlock() {this->blockname = "YEBlock";}
};

/** @class AUBlock @brief Block for AU matrix in SUSY models */
class AUBlock : public ArrayBlock<3,3> {
public:
    AUBlock() {this->blockname = "AUBlock";}
};

/** @class ADBlock @brief Block for AD matrix in SUSY models */
class ADBlock : public ArrayBlock<3,3> {
public:
    ADBlock() {this->blockname = "ADBlock";}
};

/** @class AEBlock @brief Block for AE matrix in SUSY models */
class AEBlock : public ArrayBlock<3,3> {
public:
    AEBlock() {this->blockname = "AEBlock";}
};

/** @class AlphaBlock @brief Block for alpha in SUSY models */
class AlphaBlock : public ArrayBlock<1,1> {
public:
    AlphaBlock() {this->blockname = "ALPHABlock";}
};

/** @class MSOFTBlock @brief Block for MSoft in SUSY models */
class MSOFTBlock : public MapBlock {
public:
    MSOFTBlock() {this->blockname = "MSOFTBlock";}
};
/* ------------------------------------------------------------------------------------------------
THDM BLOCKS*/

/* ------------------------------------------------------------------------------------------------
Flavor BLOCKS*/

/** @class FConstBlock @brief Block for Fconst in the flavor part */
class FMassBlock : public MapBlock {
    public:
        FMassBlock() {this->blockname = "FMASS";}
    };

/** @class FConstBlock @brief Block for Fconst in the flavor part */
class FConstBlock : public MapBlock {
public:
    FConstBlock() {this->blockname = "FCONST";}
};

/** @class FConstRatioBlock @brief Block for Fconstratio in the flavor part */
class FConstRatioBlock : public MapBlock {
public:
    FConstRatioBlock() {this->blockname = "FCONSTRATIO";}
};

/** @class FBag @brief Block for Fbag in the flavor part */
class FBag : public MapBlock {
public:
    FBag() {this->blockname = "FBAG";}
};

/** @class FLifeBlock @brief Block for Flife in the flavor part */
class FLifeBlock : public MapBlock {
public:
    FLifeBlock() {this->blockname = "FLIFE";}
};

/* ------------------------------------------------------------------------------------------------
Wilson Input BLOCKS*/

/**
 * @class WilsonBlock
 * @brief Special block for Wilson coefficients.
 * 
 * The PDG code encodes coefficient ID and QCD order.
 */
class WilsonBlock : public Block {
    // pdgCode is NNO with NN = integer ID of the coefficient as in BWilsonCoefficients enum, O is QCD order
    // pdgCode -1 is reserved to access the scale of the coefficients
    // pdgCode -2 is reserved to access the type of the coefficients
public:
    double getValue(int pdgCode) const override;

    void setValue(int pdgCode, double value, bool force = false);

    void setMode(int pdgCode, ParameterMode mode) {}

    std::map<int, double> getAllValues() override {
        return {};
    }
protected:
    // Index is QCD order
    /// Map of Wilson coefficients indexed by order
    std::map<WCoef, std::array<double, 3>> values; 
    double scale; ///< Scale of the coefficients
    int type;     ///< Type of the coefficients

};

/* ------------------------------------------------------------------------------------------------
Form Factors BLOCKS*/

/** @class BKsBlock @brief Block for BKsBlock in the form factor part */
class BKsBlock : public MapBlock {
public:
    BKsBlock() {this->blockname = "B_Ks";}
};

/** @class BllBlock @brief Block for BllBlock in the form factor part */
class BllBlock : public MapBlock {
public:
    BllBlock() {this->blockname = "B_ll";}
};

/** @class BXsBlock @brief Block for BXsBlock in the form factor part */
class BXsBlock : public MapBlock {
public:
    BXsBlock() {this->blockname = "B_Xs";}
};

/** @class BDlnuBlock @brief Block for BDlnuBlock in the form factor part */
class BDlnuBlock : public MapBlock {
public:
    BDlnuBlock() {this->blockname = "B_Dlnu";}
};

/** @class BDslnuBlock @brief Block for BDslnuBlock in the form factor part */
class BDslnuBlock : public MapBlock {
public:
    BDslnuBlock() {this->blockname = "B_Dslnu";}
};

#endif