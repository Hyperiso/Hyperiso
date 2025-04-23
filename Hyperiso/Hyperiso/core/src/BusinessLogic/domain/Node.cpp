#include "Node.h"

void ParameterNode::updateValue() {
    value = ObsParameterProxy()(lookup);
}

scalar_t ParameterNode::getValue() {
    LOG_DEBUG("ParameterNode::getValue() [", getName(), "] (", value, ")");
    return value;
}

bool ParameterNode::updateCacheFlag() {
    LOG_DEBUG("ParameterNode::updateCacheFlag() [", getName(), "]");
    if (!visited) {
        cacheValid = fpeq(value, (double)ObsParameterProxy()(lookup));
        if (!cacheValid)
            updateValue();
        visited = true;
    }
    
    return cacheValid;
}

std::string ParameterNode::getName() {
    std::stringstream ss;
    ss << "(" << lookup.block << "," << lookup.code << ")"; 
    return ss.str();
}

void ParameterNode::unvisit() {
    visited = false;
}

void ParameterNode::accept(Visitor& visitor) {
    visitor.visit(*this);
}

void OperatorNode::unvisit() {
    if (visited) {
        visited = false;
        for (auto c : children)
            c->unvisit();
    }
}

scalar_t OperatorNode::calculate() {
    updateCacheFlag();
    getValue();
    unvisit();
    return cachedValue;
}

bool OperatorNode::updateCacheFlag() {
    LOG_DEBUG("OperatorNode::updateCacheFlag() [", name, "]");
    if (!visited) {
        for (auto c : children)
            cacheValid &= c->updateCacheFlag();
        visited = true;
    }
    return cacheValid;
}

scalar_t OperatorNode::getValue() {
    LOG_DEBUG("OperatorNode::getValue() [", name, "]");
    if (!cacheValid) {
        std::vector<scalar_t> childValues;
        for (const auto& child : children) {
            childValues.push_back(child->getValue());
        }
        cachedValue = computeFunc(childValues);
        cacheValid = true;
        n_evals++;
        LOG_DEBUG("Call to OperatorNode::computeFunc [", name, "] (", cachedValue, ")");
    }
    return cachedValue;
}

void OperatorNode::addChild(const std::shared_ptr<AbstractNode>& child) {
    children.push_back(child);
}

void OperatorNode::addChildren(const std::vector<std::shared_ptr<AbstractNode>> &children) {
    for (auto c : children) {
        addChild(c);
    }
}

std::string OperatorNode::getName() {
    return name;
}

void OperatorNode::accept(Visitor& visitor) {
    visitor.visit(*this);
}