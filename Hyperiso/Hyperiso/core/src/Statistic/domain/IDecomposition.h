#ifndef IDECOMPOSITION_H
#define IDECOMPOSITION_H

#include "RNGHelper.h"

struct IDecomposition {
    virtual ~IDecomposition() = default;

    virtual Matrix factorize(const Matrix& R) = 0;
};

#endif