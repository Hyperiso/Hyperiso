#include "Observable.h"
#include "ObsEvaluator.h"
#include "Logger.h"

void Observable::evaluate() {
    this->value = ObsEvaluator::Evaluate(this);
    evaluated = this->value != complex_t(-1);
}

complex_t Observable::getValue() const {
    if (this->evaluated) {
        return this->value;
    } 
    Logger::getInstance()->error("Trying to access an unevaluated observable value.");  
    return complex_t(-1);  
}
