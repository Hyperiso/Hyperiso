#ifndef __ICONTOUREXTRACTOR_H__
#define __ICONTOUREXTRACTOR_H__

#include <vector>
#include <map>
#include <memory>

struct Contour {
    std::vector<std::pair<double, double>> points;
    std::array<double, 4> bounds;
    double level;
    bool success {false};
};


class IContourExtractor {
public:
    virtual ~IContourExtractor() = default;
    virtual void extract(std::shared_ptr<Contour>) = 0;
};


#endif // __ICONTOUREXTRACTOR_H__
