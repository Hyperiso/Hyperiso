#include "Nuisance.h"

Nuisance::Nuisance(const std::string& name, double central_value, double deviation, double sys_error)
    : name(name), central_value(central_value), deviation(deviation), sys_error(sys_error) {}

double Nuisance::get_randomized_value(int method, std::mt19937& rng) const {
    if (method == 0) return central_value;
    std::normal_distribution<> dist(central_value, deviation);
    return dist(rng);
}