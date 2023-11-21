#include <complex>
#include <string>

typedef std::complex<double> complex_t;

enum class RenormalizationScheme {
    POLE, MSBAR, DRBAR, ONE_S, KIN, INVARIANT, MOM, SMOM
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

// TODO: Implement constructor forwarding
template <typename T>
class SchemeDependentElement : public TypedElement<T> {
private:
    const RenormalizationScheme scheme;

public:
    SchemeDependentElement(const std::string& id, const T& value, RenormalizationScheme scheme) : TypedElement<T>(id, value), scheme(scheme) {}
    RenormalizationScheme getScheme() const { return this->scheme; }
};

template <typename T>
class ScaleDependentElement : public TypedElement<T> {
private:
    double Q;

public:
    ScaleDependentElement(const std::string& id, const T& value, double scale) : TypedElement<T>(id, value), Q(scale) {}
    double getScale() const { return this->Q; }
};

template <typename T>
class GeneralElement : public SchemeDependentElement<T>, public ScaleDependentElement<T> {

public:
    GeneralElement(const std::string& id, const T& value, RenormalizationScheme scheme, double scale) 
        : SchemeDependentElement<T>(id, value, scheme), 
        ScaleDependentElement<T>(id, value, scale) {}
};