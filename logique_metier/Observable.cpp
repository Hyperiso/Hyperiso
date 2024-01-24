#include "Observable.h"
#include "ObsEvaluator.h"

void Observable::evaluate() {
    this->value = ObsEvaluator::Evaluate(this);
    if (this->value != complex_t(-1)) {
        evaluated = true;
    } 
}

complex_t Observable::getValue() const {
    if (this->evaluated) {
        return this->value;
    } 
    // Log an error    
}