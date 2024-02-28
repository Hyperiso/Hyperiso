#include "Observable.h"
#include "ObsEvaluator.h"
#include "./Core/Logger.cpp"

void Observable::evaluate() {
    this->value = ObsEvaluator::Evaluate(this);
    evaluated = this->value != complex_t(-1);
}

complex_t Observable::getValue() const {
    if (this->evaluated) {
        return this->value;
    } 
    auto log = Logger::getInstance();    
    log.error("Trying to access an unevaluated observable value.");
}
