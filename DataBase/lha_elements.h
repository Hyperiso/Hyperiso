#ifndef HYPERISO_LHA_ELEMENTS_H
#define HYPERISO_LHA_ELEMENTS_H

#include <complex>
#include <string>

typedef std::complex<double> complex_t;

enum class RenormalizationScheme {
    POLE, MSBAR, DRBAR, ONE_S, KIN, INVARIANT, MOM, SMOM, NONE
};

class AbstractElement {
protected:
    const std::string id;

public:
    AbstractElement(const std::string& id) : id(id) {}
    std::string getId() const { return this->id; }
    virtual std::string toString() const = 0;
    virtual ~AbstractElement() = default;
};

template <typename T>
class TypedElement : public AbstractElement {
private:
    T value;

public:
    TypedElement(const std::string& id, const T& value) : AbstractElement(id), value(value) {}
    T getValue() const { return this->value; }

    std::string toString() const override {
        std::stringstream stream;
        stream << this->getId() << '\t' << this->getValue() << "\n";
        return stream.str();
    }
};

template <typename T>
class GeneralElement : public TypedElement<T> {
private:
    const RenormalizationScheme scheme;
    double Q;

public:
    GeneralElement(const std::string& id, const T& value, RenormalizationScheme scheme, double scale) 
        : TypedElement<T>(id, value), scheme(scheme), Q(scale) {}
    RenormalizationScheme getScheme() const { return this->scheme; }
    double getScale() const { return this->Q; }
};

template <typename T>
class SchemeDependentElement : public GeneralElement<T> {

public:
    SchemeDependentElement(const std::string& id, const T& value, RenormalizationScheme scheme) 
        : GeneralElement<T>(id, value, scheme, 0) {}
};

template <typename T>
class ScaleDependentElement : public GeneralElement<T> {

public:
    ScaleDependentElement(const std::string& id, const T& value, double scale) 
        : GeneralElement<T>(id, value, RenormalizationScheme::NONE, scale) {}
};

#endif // HYPERISO_LHA_ELEMENTS_H