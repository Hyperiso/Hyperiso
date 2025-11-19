#ifndef IDECOMPOSITION_H
#define IDECOMPOSITION_H

#include "RNGHelper.h"
#include "Include.h"
#include <map>

struct IDecomposition {
    virtual ~IDecomposition() = default;

    virtual Matrix factorize(const Matrix& R) = 0;
    virtual std::map<ParamId, std::map<ParamId, double>> factorize(const std::map<ParamId, std::map<ParamId, double>>& R) = 0;
};

#endif