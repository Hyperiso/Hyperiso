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
// #include "JsonParameters.h"
#include "Parameter.h"

/**
 * @class Block
 * @brief Abstract base class for parameter blocks.
 */
class Block {
public:
    /**
     * @brief Retrieves the value associated with a given PDG code (int).
     * @param id The PDG code of the parameter.
     * @return The parameter value.
     */
    virtual double getValue(LhaID id) const = 0;

    /**
     * @brief Retrieves the value associated with a given PDG code (int).
     * @param id The PDG code of the parameter.
     * @return The parameter.
     */
    virtual Parameter getParameter(LhaID id) const = 0;

    /**
     * @brief Sets the value of a parameter.
     * @param id The PDG code of the parameter.
     * @param value The new value to set.
     * @param force If true, forces the update.
     */
    virtual void setValue(LhaID id, double value, bool force = false) = 0;

    /**
     * @brief Sets the value of a parameter.
     * @param id The PDG code of the parameter.
     * @param source The source parameter to set.
     */
    virtual void setParameter(LhaID id, const Parameter& source) = 0;

    /**
     * @brief Sets the value of a parameter.
     * @param id The PDG code of the parameter.
     * @param std_stat The new deviation to set.
     * @param std_syst The new deviation to set.
     * @param force If true, forces the update.
     */
    virtual void setDeviation(LhaID id, double std_stat, double std_syst, bool force = false) = 0;

    /**
     * @brief Sets the mode of a parameter.
     * @param id The PDG code of the parameter.
     * @param mode The mode to set.
     */
    virtual void setMode(LhaID id, ParameterMode mode) = 0;

    /**
     * @brief Retrieves all parameter values.
     * @return A map of PDG codes to parameter values.
     */
    virtual std::map<LhaID, double> getAllValues() = 0;

    /**
     * @brief Retrieves all parameter ids.
     * @return A vector of LhaIDs of stored parameters.
     */
    virtual std::vector<LhaID> getAllIDs() = 0;

    /**
     * @brief Retrieves all parameters with their ids.
     * @return A map of LhaIDs and Parameters.
     */
    virtual const std::map<LhaID, Parameter>& getItems() = 0;

    /**
     * @brief Retrieves all parameter ids.
     * @return A vector of LhaIDs of stored parameters.
     */
    virtual bool hasID(LhaID id) = 0;

    /**
     * @brief Removes parameter from block.
     * @param id Id of the parameter to remove.
     */
    virtual void remove_parameter(LhaID id) = 0;

    Block() = default;

    Block(std::shared_ptr<Block> other) { this->copy(other); };

    void addObserver(std::shared_ptr<Block> observer);

    void notifyObservers();

    virtual void copy(std::shared_ptr<Block> other) = 0;

    virtual void update() = 0;
    /**
     * @brief Virtual destructor.
     */
    virtual ~Block() = default;

    /// Name of the block
    std::string blockname{};
protected:
    std::vector<std::shared_ptr<Block>> observers;
};

/**
 * @class MapBlock
 * @brief A class for  storing blocks that have only 1 id.
 */
class MapBlock : public Block {
public:

    MapBlock() = default;
    MapBlock(std::shared_ptr<Block> other);
    
    double getValue(LhaID id) const override;
    void setValue(LhaID id, double value, bool force = false) override;
    void setDeviation(LhaID id, double std_stat, double std_syst, bool force = false) override;
    void setMode(LhaID id, ParameterMode mode) override;
    std::map<LhaID, double> getAllValues() override;
    std::vector<LhaID> getAllIDs() override;
    bool hasID(LhaID id) override;
    void copy(std::shared_ptr<Block> other) override;
    const std::map<LhaID, Parameter>& getItems() override { return this->values; };
    Parameter getParameter(LhaID id) const override;
    void setParameter(LhaID id, const Parameter& source) override;
    void remove_parameter(LhaID id) override;

    void update() override {
        recalculate();
    }

    ~MapBlock() { notifyObservers(); }

protected:
    /// Map of PDG codes to parameters
    std::map<LhaID, Parameter> values;

    void recalculate() {
        std::cout << "Error map block cannot do that" << std::endl;
    }

};

class DependentBlock : public MapBlock, public std::enable_shared_from_this<DependentBlock> {
public:
    explicit DependentBlock(std::shared_ptr<Block> source) 
        : sourceBlock(source) {}

    void init();

    void update();

protected:
    std::shared_ptr<Block> sourceBlock;

    virtual void recalculate();
};

class GaugeBlock: public DependentBlock {
protected:
    void recalculate() override;
};

class ReCKMBlock: public DependentBlock {
protected:
    void recalculate() override;
};

class ImCKMBlock: public DependentBlock {
protected:
    void recalculate() override;
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
    double getValue(LhaID pdgCode) const override;

    void setValue(LhaID pdgCode, double value, bool force = false);

    void setMode(LhaID pdgCode, ParameterMode mode) {}

    std::map<LhaID, double> getAllValues() override {
        return {};
    }

protected:
    // Index is QCD order
    /// Map of Wilson coefficients indexed by order
    std::map<WCoef, std::array<double, 3>> values; 
    double scale; ///< Scale of the coefficients
    int type;     ///< Type of the coefficients

};

#endif