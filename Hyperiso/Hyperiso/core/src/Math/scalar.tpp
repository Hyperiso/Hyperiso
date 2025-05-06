#include "scalar.h"

template<typename T>
requires std::convertible_to<T, std::complex<double>>
scalar_t& scalar_t::operator+=(const T& rhs) {
    // *this = static_cast<complex_t>(*this) + static_cast<complex_t>(rhs);
    complex_t::operator+=(static_cast<const complex_t&>(rhs));
    return *this;
}

template<typename T>
requires std::convertible_to<T, std::complex<double>>
scalar_t operator+(scalar_t lhs, const T& rhs) {
    lhs += rhs;
    return lhs;
}

template<typename T>
requires std::convertible_to<T, std::complex<double>>
scalar_t operator+(T lhs, const scalar_t& rhs) {
    lhs += rhs;
    return lhs;
}

template<typename T>
requires std::convertible_to<T, std::complex<double>>
scalar_t& scalar_t::operator-=(const T& rhs) {
    *this = static_cast<complex_t>(*this) - static_cast<complex_t>(rhs);
    return *this;
}

template<typename T>
requires std::convertible_to<T, std::complex<double>>
scalar_t operator-(scalar_t lhs, const T& rhs) {
    lhs -= rhs;
    return lhs;
}

template<typename T>
requires std::convertible_to<T, std::complex<double>>
scalar_t operator-(T lhs, const scalar_t& rhs) {
    return scalar_t(static_cast<std::complex<double>>(lhs) - static_cast<std::complex<double>>(rhs));

}

template<typename T>
requires std::convertible_to<T, std::complex<double>>
scalar_t& scalar_t::operator*=(const T& rhs)
{
    *this = static_cast<complex_t>(*this) * static_cast<complex_t>(rhs);
    return *this;
}

template<typename T>
requires std::convertible_to<T, std::complex<double>>
scalar_t operator*(scalar_t lhs, const T& rhs) {
    lhs *= rhs;
    return lhs;
}

template<typename T>
requires std::convertible_to<T, std::complex<double>>
scalar_t operator*(T lhs, const scalar_t& rhs) {
    lhs *= rhs;
    return lhs;
}

template<typename T>
requires std::convertible_to<T, std::complex<double>>
scalar_t& scalar_t::operator/=(const T& rhs) {
    *this = static_cast<complex_t>(*this) / static_cast<complex_t>(rhs);
    return *this;
}


template<typename T>
requires std::convertible_to<T, std::complex<double>>
scalar_t operator/(scalar_t lhs, const T& rhs) {
    lhs /= rhs;
    return lhs;
}

template<typename T>
requires std::convertible_to<T, std::complex<double>>
scalar_t operator/(T lhs, const scalar_t& rhs) {
    lhs /= rhs;
    return lhs;
}

template <typename T>
requires std::is_integral_v<T>
scalar_t pow(const scalar_t& base, T exp) {
    return pow(base, static_cast<double>(exp));
}