#include <string>
#include <cmath>
#include <unordered_map>
#include "Logger.h"

enum DistributionType {
    Gaussian = 1,
    Flat = 2
};

struct indnuis {
    double cent, dev;
    DistributionType type;
    std::string name;

    indnuis(double c, double d, DistributionType t, std::string n)
        : cent(c), dev(d), type(t), name(std::move(n)) {}
};

class Nuisance {
    
    std::unordered_map<std::string, indnuis> parameters;

public:
    indnuis alphas_MZ, mass_b, mass_c, mass_s, mass_top_pole, mass_h0;
    Nuisance()
        : alphas_MZ(0.1181, 0.0011, Gaussian, "alphas_MZ"),
          mass_b(4.18, 0.04, Gaussian, "mass_b"),
          mass_c(1.27, 0.003, Gaussian, "mass_c"),
          mass_s(0.096, 0.008, Gaussian, "mass_s"),
          mass_top_pole(173.34, sqrt(0.27*0.27+0.71*0.71), Gaussian, "mass_top_pole"),
          mass_h0(125.09, 0.24, Gaussian, "mass_h0") {
          }

    void setParameter(const std::string& key, const indnuis& value) {
        parameters[key] = value;
    }

    const indnuis& getParameter(const std::string& key) const {
        auto it = parameters.find(key);
        if (it != parameters.end()) {
            return it->second;
        }
        LOG_ERROR("Parameter not found");
    }

    bool hasParameter(const std::string& key) const {
        return parameters.find(key) != parameters.end();
    }
};


class NuisanceBuilder {
private:
    Nuisance n;

public:
    NuisanceBuilder& setAlphasMZ(double cent, double dev) {
        n.alphas_MZ = indnuis(cent, dev, Gaussian, "alphas_MZ");
        return *this;
    }

    NuisanceBuilder& setMassB(double cent, double dev) {
        n.mass_b = indnuis(cent, dev, Gaussian, "mass_b");
        return *this;
    }


    Nuisance build() {
        return n;
    }
};

