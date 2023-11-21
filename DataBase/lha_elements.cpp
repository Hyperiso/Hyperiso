#include <complex>
#include <string>

typedef std::complex<double> complex_t;

// Interface de base pour les données individueles, template pour avoir une value réelle (double) ou complexe (std::complex<double>)
class LhaElement {
private:
    const int id;
    complex_t value;

public:
    LhaElement(int id, complex_t value) : id(id), value(value) {}
    int getId() const { return this->id; }
    complex_t getValue() const { return this->value; }
};

class LhaReal : public LhaElement {
public:
    LhaReal(int id, double value) : LhaElement(id, static_cast<complex_t>(value)) {}
};