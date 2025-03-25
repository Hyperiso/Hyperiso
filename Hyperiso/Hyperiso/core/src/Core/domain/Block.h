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

class Block;
class DependentBlock;
typedef std::function<void(const std::unordered_map<std::string, std::shared_ptr<Block>>&, std::shared_ptr<DependentBlock>)> DepUpdateFunc;

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
    void store(const LhaID& id, Parameter&& param) override;
    void assign(const LhaID& key, Parameter&& param) override;
    void store_or_assign(const LhaID& key, Parameter&& param) override;
    bool contains(const LhaID& key) const override;
    Parameter& retrieve(const LhaID& id) override;
    void remove(const LhaID& key) override;

    std::unordered_set<LhaID> getAllIDs();
    const std::map<LhaID, Parameter>& getItems() { return this->items; };
    void set_owner(ParameterType type);

    void addObserver(std::shared_ptr<Block> observer);
    void removeObserver(std::shared_ptr<Block> observer);
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
    explicit DependentBlock(std::unordered_map<std::string, std::shared_ptr<Block>> sources, DepUpdateFunc recalculateFunc) 
        : sourceBlocks(std::move(sources)), recalculateLambda(std::move(recalculateFunc)) {}

    bool dependsOn(const std::string& blockName) {
        return sourceBlocks.contains(blockName);
    }

    void init() {
        self = shared_from_this();
        if (self) {
            for (auto src : sourceBlocks){
                src.second->addObserver(self);   
            }
        } else {
            std::cerr << "Error: DependentBlock must be created with std::make_shared!" << std::endl;
        }
    }

    void update() override {
        if (recalculateLambda 
            && std::all_of(sourceBlocks.begin(), sourceBlocks.end(), 
                           [](std::pair<std::string, std::shared_ptr<Block>> block) { return block.second; })) 
        {
            if (auto self = shared_from_this()) { 
                recalculateLambda(sourceBlocks, self);
            } else {
                std::cerr << "Error: shared_from_this() failed in update()" << std::endl;
            }
        }
    }

    ~DependentBlock() {
        LOG_INFO("Destruct dependentBlock at", self.get());
        if (self) {
            for (auto src : sourceBlocks){
                src.second->removeObserver(self);   
            }
        }
    }

private:
    std::shared_ptr<DependentBlock> self;
    std::unordered_map<std::string, std::shared_ptr<Block>> sourceBlocks;
    DepUpdateFunc recalculateLambda;
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