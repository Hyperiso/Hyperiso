#ifndef __ICONTOUREXTRACTOR_H__
#define __ICONTOUREXTRACTOR_H__

#include "Math.h"
#include <vector>
#include <set>
#include <memory>
#include <functional>

struct ContourRequest {
    double level;
    std::array<double, 4> bounds;
    std::array<fit_app::ParameterDefinition, 2> p_defs;
    std::size_t resolution = 40;
};

struct Contour {
    std::set<Path> paths;
    double level;
    bool success = false;
};

class IContourExtractor {
public:
    virtual ~IContourExtractor() = default;
    virtual Contour extract(const ScalarField2D& field, const ContourRequest& cr) = 0;
};


#endif // __ICONTOUREXTRACTOR_H__
