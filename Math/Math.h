double Li2(double x);
double Cl2(double x);
double H2(double x, double y);
double B(double m1, double m2, double Q);
template <typename T> int sgn(T val) {
    return (T(0) < val) -(val<T(0));
}