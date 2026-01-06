#ifndef IOBS_CORRELATION_PROXY_H
#define IOBS_CORRELATION_PROXY_H

#include <concepts>

#include "Math.h"
#include "ParameterProvider.h"

template<typename T, typename V, typename U, typename W>
class IObsCorrelationProxy {
public:
    virtual double operator()(const T& pid_1, const T& pid_2, W type) = 0;

    virtual double operator()(const V& pid_1, const V& pid_2, W type) = 0;

    virtual double operator()(const U& pid_1, const U& pid_2, W type) = 0;
};

#endif