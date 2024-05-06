#include <cmath>

double I0(double x) {
    double y;
    double absx = std::fabs(x);

    if (absx < 3.75) {
        y = x / 3.75;
        y = y * y;
        return 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.0360768 + y * 0.0045813)))));
    } else {
        y = 3.75 / absx;
        return (std::exp(absx) / std::sqrt(absx)) * (0.39894228 + y * (0.01328592 + y * (0.00225319 + y * (-0.00157565 + y * (0.00916281 + y * (-0.02057706 + y * (0.02635537 + y * (-0.01647633 + y * 0.00392377))))))));
    }
}

// J. Olivares et al 2018 J. Phys.: Conf. Ser. 1043 012003
// double I0(double x) {
//     double xs = x * x;
//     return std::cosh(x) * (1 + 0.24273 * xs) / (std::pow(1 + 0.25 * xs, 0.25) * (1 + 0.43023 * xs));
// }


double I1(double x) {
    double y, tmp;
    double I1 = 0.0;
    double absx = std::fabs(x);

    if (absx < 3.75) {
        y = x / 3.75;
        y = y * y;
        I1 = absx * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934 + y * (0.02658733 + y * (0.00301532 + y * 0.00032411))))));
    } else {
        y = 3.75 / absx;
        tmp = 0.02282967 + y * (-0.02895312 + y * (0.01787654 - y * 0.00420059));
        tmp = 0.39894228 + y * (-0.03988024 + y * (-0.00362018 + y * (0.00163801 + y * (-0.01031555 + y * tmp))));
        I1 *= (std::exp(absx) / std::sqrt(absx));
    }

    if (x < 0.0) {
        return -I1;
    } else {
        return I1;
    }
}

double K0(double x) {
    double y, result;

    if (x <= 2.0) {
        y = x * x / 4.0;
        result = (-std::log(x / 2.0) * I0(x)) + (-0.57721566 + y * (0.42278420 + y * (0.23069756 + y * (0.03488590e-1 + y * (0.00262698e-2 + y * (0.00010750e-3 + y * 0.74e-5))))));
    } else {
        y = 2.0 / x;
        result = (std::exp(-x) / std::sqrt(x)) * (1.25331414 + y * (-0.07832358e-1 + y * (0.02189568e-1 + y * (-0.01062446e-1 + y * (0.00587872e-2 + y * (-0.00251540e-2 + y * 0.00053208e-3))))));
    }

    return result;
}


double K1(double x) {
    double y, result;

    if (x <= 2.0) {
        y = x * x / 4.0;
        result = (std::log(x / 2.0) * I1(x)) + (1.0 / x) * (1.0 + y * (0.15443144 + y * (-0.67278579 + y * (-0.18156897 + y * (-0.01919402e-1 + y * (-0.00110404e-2 + y * (-0.4686e-4)))))));
    } else {
        y = 2.0 / x;
        result = (std::exp(-x) / std::sqrt(x)) * (1.25331414 + y * (0.23498619 + y * (-0.03655620e-1 + y * (0.01504268e-1 + y * (-0.00780353e-2 + y * (0.00325614e-2 + y * (-0.0068245e-3)))))));
    }

    return result;
}


double K2(double x)
/* calculates the Modified Bessel functions of second type and of order 2 of x */
{
	return K0(x)+2./x*K1(x);
}


double K3(double x)
/* calculates the Modified Bessel functions of second type and of order 3 of x */
{
	return K1(x)+4./x*K2(x);
}


double K4(double x)
/* calculates the Modified Bessel functions of second type and of order 4 of x */
{
	return K2(x)+6./x*K3(x);
}


double Lbessel(double x)
/* calculates L(x)=K2(x)/x */
{
	return K2(x)/x;
}


double Mbessel(double x)
/* calculates M(x)=(3*K3(x)+K1(x))/4x */
{
	return (0.75*K3(x)+0.25*K1(x))/x;
}


double Nbessel(double x)
/* calculates N(x)=(K4(x)+K2(x))/2x */
{
	return (0.5*K4(x)+0.5*K2(x))/x;
}



double K0exp(double x, double z) {
  /* calculates the extended Modified Bessel functions of second type and of order 0 of x and z */
    double y, result;

    if (x <= 2.0) {
        y = x * x / 4.0;
        result = ((-std::log(x / 2.0) * I0(x)) + (-0.57721566 + y * (0.42278420 + y * (0.23069756 + y * (0.03488590e-1 + y * (0.0262698e-2 + y * (0.0010750e-3 + y * 0.74e-5))))))) * std::exp(z) * std::sqrt(z);
    } else {
        y = 2.0 / x;
        result = (std::exp(z - x) / std::sqrt(x / z)) * (1.25331414 + y * (-0.07832358e-1 + y * (0.02189568e-1 + y * (-0.01062446e-1 + y * (0.00587872e-2 + y * (-0.00251540e-2 + y * 0.53208e-3))))));
    }

    return result;
}


double K1exp(double x, double z) {
    double y, result;

    if (x <= 2.0) {
        y = x * x / 4.0;
        result = ((std::log(x / 2.0) * I1(x)) + (1.0 / x) * (1.0 + y * (0.15443144 + y * (-0.67278579 + y * (-0.18156897 + y * (-0.01919402e-1 + y * (-0.00110404e-2 + y * -0.4686e-4))))))) * std::exp(z) * std::sqrt(z);
    } else {
        y = 2.0 / x;
        result = (std::exp(z - x) / std::sqrt(x / z)) * (1.25331414 + y * (0.23498619 + y * (-0.03655620e-1 + y * (0.01504268e-1 + y * (-0.00780353e-2 + y * (0.00325614e-2 + y * -0.68245e-3))))));
    }

    return result;
}

double K2exp(double x,double z)
/* calculates the extended Modified Bessel functions of second type and of order 2 of x and z */
{
	return K0exp(x,z)+2./x*K1exp(x,z);
}






