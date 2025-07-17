#ifndef OBS_NODE_H
#define OBS_NODE_H

#include <iostream>
#include <vector>
#include <memory>
#include <functional>
#include <unordered_set>
#include <map>
#include <chrono>

#include "ObsParameterProxy.h"
#include "General.h"
#include "Math.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration;

// Node interface, purely virtual
class AbstractNode {

public:
    virtual ~AbstractNode() = default;

    // Access the node's value
    virtual scalar_t getValue() = 0;

    // Make cache invalid and propagate the information upwards 
    virtual bool updateCacheFlag() = 0;

    virtual void unvisit() = 0;

    // Accept a visitor (Visitor Pattern)
    virtual void accept(class Visitor& visitor) = 0;

    virtual std::string getName() = 0;
};

// Parameter node (leaf)
class ParameterNode : public AbstractNode {
private:
    ParamId lookup;                // Parameter id to cache
    scalar_t value;                  // Cached value
    bool cacheValid;
    bool visited;

public:
    inline ParameterNode(ParamId param_id) : lookup(param_id), value(ObsParameterProxy()(param_id)), cacheValid(true), visited(false) {}

    // Updates cached value by reading the parameter structure
    void updateValue();

    // Always return the up-to-date value.
    scalar_t getValue() override;

    // Checks whether the cached value is still the same as the one stored in the parameters
    bool updateCacheFlag() override;

    void unvisit() override;

    std::string getName() override;

    void accept(Visitor& visitor) override;
};

// === Classe pour les opérateurs (Composite Pattern) ===
class OperatorNode : public AbstractNode {
private:
    std::string name;
    std::function<scalar_t([[maybe_unused]] const std::vector<scalar_t>&)> computeFunc; // Fonction de calcul
    std::vector<std::shared_ptr<AbstractNode>> children;           // Enfants
    scalar_t cachedValue;
    bool cacheValid;
    bool visited;
    size_t n_evals;

public:
    template <typename Callable>
    OperatorNode(std::string name, Callable&& func)
        : name(std::move(name)), computeFunc(std::forward<Callable>(func)), cachedValue(0.0), cacheValid(false), visited(false), n_evals(0) {}

    void addChild(const std::shared_ptr<AbstractNode>& child);
    void addChildren(const std::vector<std::shared_ptr<AbstractNode>>& children);

    scalar_t getValue() override;

    // Cache remains valid iff all children caches do
    bool updateCacheFlag() override;

    // Recursively unvisit all children nodes
    void unvisit() override;

    // Complete calculation including cache checking
    scalar_t calculate();

    std::string getName();

    size_t get_n_evals() { return n_evals; }

    void accept(Visitor& visitor) override;
};

// === Visitor Pattern ===
class Visitor {
public:
    virtual void visit(ParameterNode& node) = 0;
    virtual void visit(OperatorNode& node) = 0;
};

// === Exemple de visiteur : affichage des noeuds ===
class PrintVisitor : public Visitor {
public:
    void visit(ParameterNode& node) override {
        std::cout << "mmh" << std::endl;
        std::cout << "ParameterNode " << node.getName() << " : " << node.getValue() << "\n";
    }

    void visit(OperatorNode& node) override {
        scalar_t value = node.calculate();
        std::cout << "OperatorNode " << node.getName() << " : " << value << "\n";
    }
};

#endif