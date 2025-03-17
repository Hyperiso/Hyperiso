/**
 * @file Blocks.h
 * @brief Defines various block classes for parameter management.
 * 
 * This file declares a hierarchy of block classes used to store and manage parameter values in the differents 
 * Parameters instances.
 */
#if !defined(BLOCK_H)
#define BLOCK_H

#include <array>
#include <functional>
#include "Include.h"
#include "Parameter.h"
#include "IStorage.h"

/**
 * @class Block
 * @brief A class for storing parameter blocks.
 */
class Block : public IStorage<LhaID, Parameter> {
public:
    std::string blockname {""};

    Block() = default;
    Block(std::shared_ptr<Block> other);

    // Interface methods
    const Parameter& retrieve(const LhaID& id) const override;
    void store(const LhaID& id, Parameter&& param) override;
    void remove(const LhaID& key) override;
    bool contains(const LhaID& key) const override;
    void update(const LhaID& key, Parameter&& param) override;

    std::unordered_set<LhaID> getAllIDs();
    const std::map<LhaID, Parameter>& getItems() { return this->items; };

    void addObserver(std::shared_ptr<Block> observer);
    void notifyObservers();
    virtual void update() {}

    void copy(std::shared_ptr<Block> other);

    ~Block() { notifyObservers(); }

protected:
    std::vector<std::shared_ptr<Block>> observers;
    std::map<LhaID, Parameter> items;
};

class DependentBlock : public Block, public std::enable_shared_from_this<DependentBlock> {
public:
    explicit DependentBlock(std::shared_ptr<Block> source, std::function<void(std::shared_ptr<Block>, std::shared_ptr<DependentBlock>)> recalculateFunc) 
        : sourceBlock(source), recalculateLambda(std::move(recalculateFunc)) {}

    void init() {
        if (auto self = shared_from_this()) {
            sourceBlock->addObserver(self);
        } else {
            std::cerr << "Error: DependentBlock must be created with std::make_shared!" << std::endl;
        }
    }

    void update() override {
        if (recalculateLambda && sourceBlock) {
            if (auto self = shared_from_this()) { 
                recalculateLambda(sourceBlock, self);
            } else {
                std::cerr << "Error: shared_from_this() failed in update()" << std::endl;
            }
        }
    }

private:
    std::shared_ptr<Block> sourceBlock;
    std::function<void(std::shared_ptr<Block>, std::shared_ptr<DependentBlock>)> recalculateLambda;
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
    double getValue(LhaID pdgCode) const;
    void setValue(LhaID pdgCode, double value, bool force = false);

protected:
    // Index is QCD order
    /// Map of Wilson coefficients indexed by order
    std::map<WCoef, std::array<double, 3>> values; 
    double scale; ///< Scale of the coefficients
    int type;     ///< Type of the coefficients

};

#endif