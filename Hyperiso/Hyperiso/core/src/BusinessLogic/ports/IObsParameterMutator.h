#ifndef IOBS_PARAMETER_MUTATOR_H
#define IOBS_PARAMETER_MUTATOR_H

#include "Include.h"
#include "Math.h"

template <typename T, typename U, typename V>
class IObsParameterMutator {
public:
    virtual ~IObsParameterMutator() = default;
    
    virtual void set(const T&, U) = 0;
    virtual void shift(const T&, U) = 0;
    virtual void change_mode(const T&, V) = 0;
};

#endif // __IPARAMMODIFIER_H__
