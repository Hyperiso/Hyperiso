#include "Node.h"

void ParameterNode::updateValue() {
    value = Parameters::Get(lookup);
}

scalar_t ParameterNode::getValue() {
    // std::cout << "ParameterNode::getValue() [" << getName() << "] (" << value << ")" << std::endl;
    return value;
}

bool ParameterNode::updateCacheFlag() {
    // std::cout << "ParameterNode::updateCacheFlag() [" << getName() << "] ";
    if (!visited) {
        cacheValid = fpeq(value, Parameters::Get(lookup));
        if (!cacheValid)
            updateValue();
        visited = true;
        // std::cout << (cacheValid ? "(o)" : "(!)") << std::endl;
    } else { /* std::cout << "(x)" << std::endl; */ }
    
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
    // std::cout << "OperatorNode::updateCacheFlag() [" << name << "]" << std::endl;
    if (!visited) {
        for (auto c : children)
            cacheValid &= c->updateCacheFlag();
        visited = true;
    }
    return cacheValid;
}

scalar_t OperatorNode::getValue() {
    // std::cout << "OperatorNode::getValue() [" << name << "]" << std::endl;
    if (!cacheValid) {
        std::vector<scalar_t> childValues;
        for (const auto& child : children) {
            childValues.push_back(child->getValue());
        }
        cachedValue = computeFunc(childValues);
        cacheValid = true;
        std::cout << "Call to OperatorNode::computeFunc [" << name << "] (" << cachedValue << ")" << std::endl;
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