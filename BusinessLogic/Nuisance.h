#ifndef NUISANCE_H
#define NUISANCE_H

#include <string>
#include <random>

class Nuisance {
public:
    std::string name;
    double central_value;
    double deviation;
    double sys_error;

    Nuisance(const std::string& name, double central_value, double deviation, double sys_error);
    Nuisance() = default;
    double get_randomized_value(int method, std::mt19937& rng) const;
};

#endif // NUISANCE_H