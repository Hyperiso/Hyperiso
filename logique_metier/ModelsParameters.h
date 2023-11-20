#ifndef MODEL_PARAMETERS_H
#define MODEL_PARAMETERS_H

class ModelParameters {
public:
    ModelParameters(double mass, double coupling) : particleMass(mass), couplingConstant(coupling) {}

    double getParticleMass() const { return particleMass; }
    double getCouplingConstant() const { return couplingConstant; }

private:
    double particleMass;
    double couplingConstant;
};

#endif // MODEL_PARAMETERS_H
