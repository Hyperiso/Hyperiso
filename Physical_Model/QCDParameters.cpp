#include "./QCDParameters.h"
#include <iostream>

QCDParameters::QCDParameters() {
    mass_Z = 91.1699982;
    alphas_MZ = 0.117200002;
    
}

double QCDParameters::alphasRunning(double Q, double Lambda, int nf) {
    double beta0 = 11.-2./3.*nf;
    double beta1=51.-19./3.*nf;
    double beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

    return 4.*pi/beta0/log(pow(Q/Lambda,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda,2.)))/log(pow(Q/Lambda,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda,2.)),2.)*(pow(log(log(pow(Q/Lambda,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
}

double QCDParameters::DichotomieLambda(double alpha_running, double Q, int nf) {
    
    int _ = nf;
    nf = 5;
    double beta0 = 11.-2./3.*nf;
    double beta1=51.-19./3.*nf;
    double beta2=2857.-5033.*nf/9.+325./27.*nf*nf;
    Lambda5 = 0.226;
    
    
    double alphas_running=4.*pi/beta0/log(pow(Q/Lambda5,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda5,2.)))/log(pow(Q/Lambda5,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda5,2.)),2.)*(pow(log(log(pow(Q/Lambda5,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
    nf = _;
    beta0 = 11.-2./3.*nf;
    beta1=51.-19./3.*nf;
    beta2=2857.-5033.*nf/9.+325./27.*nf*nf;
    std::cout << alphas_running << std::endl;
    double Lambda_min=1.e-3;
    double Lambda_max=1.;
    double alphas_min=0.;
    double alphas_max=0;
    double alphas_moy = 0;
    double Lambda_moy = 0;
    while((fabs(1.-alphas_min/alphas_running)>=1.e-4)&&(fabs(1.-Lambda_min/Lambda_max)>1.e-5))
    {
        alphas_min=4.*pi/beta0/log(pow(mtop/Lambda_min,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda_min,2.)))/log(pow(mtop/Lambda_min,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda_min,2.)),2.)*(pow(log(log(pow(mtop/Lambda_min,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

        alphas_max=4.*pi/beta0/log(pow(mtop/Lambda_max,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda_max,2.)))/log(pow(mtop/Lambda_max,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda_max,2.)),2.)*(pow(log(log(pow(mtop/Lambda_max,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

        Lambda_moy=(Lambda_min+Lambda_max)/2.;
        alphas_moy=4.*pi/beta0/log(pow(mtop/Lambda_moy,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(mtop/Lambda_moy,2.)))/log(pow(mtop/Lambda_moy,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(mtop/Lambda_moy,2.)),2.)*(pow(log(log(pow(mtop/Lambda_moy,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));

        if((alphas_running>=alphas_min)&&(alphas_running<=alphas_moy))
            Lambda_max=Lambda_moy;
        else Lambda_min=Lambda_moy;
    }

    Lambda6=Lambda_min;
    
    if(fabs(1.-Lambda_min/Lambda_max)<=1.e-5)
    {
        Lambda5=-1.;
        return -1.;
    }
    return Lambda_min;
}


// double QCDParameters::calculateLambda(double MZ, double alphas_MZ, int nf) {
//         const double pi = 3.14159265358979323846;
//         double beta0 = 11.0 - (2.0 / 3.0) * nf;
//         double beta1 = 51.0 - (19.0 / 3.0) * nf;
//         double beta2 = 2857.0 - (5033.0 / 9.0) * nf + (325.0 / 27.0) * nf * nf;

//         double Lambda_min = 1e-3;
//         double Lambda_max = 1.0;
//         double alphas_min, alphas_moy;

//         while (true) {
//             double log_MZ_Lambda_min = log(pow(MZ / Lambda_min, 2));
//             alphas_min = 4.0 * pi / (beta0 * log_MZ_Lambda_min) *
//                          (1.0 - 2.0 * beta1 / (beta0 * beta0) * log(log_MZ_Lambda_min) / log_MZ_Lambda_min +
//                           4.0 * beta1 * beta1 / (beta0 * beta0 * beta0 * beta0 * log_MZ_Lambda_min * log_MZ_Lambda_min) *
//                           (pow(log(log_MZ_Lambda_min) - 0.5, 2) + beta2 * beta0 / (8.0 * beta1 * beta1) - 5.0 / 4.0));

//             Lambda_max = (Lambda_min + Lambda_max) / 2.0;
//             double log_MZ_Lambda_max = log(pow(MZ / Lambda_max, 2));
//             alphas_moy = 4.0 * pi / (beta0 * log_MZ_Lambda_max) *
//                          (1.0 - 2.0 * beta1 / (beta0 * beta0) * log(log_MZ_Lambda_max) / log_MZ_Lambda_max +
//                           4.0 * beta1 * beta1 / (beta0 * beta0 * beta0 * beta0 * log_MZ_Lambda_max * log_MZ_Lambda_max) *
//                           (pow(log(log_MZ_Lambda_max) - 0.5, 2) + beta2 * beta0 / (8.0 * beta1 * beta1) - 5.0 / 4.0));

//             if ((alphas_MZ >= alphas_min) && (alphas_MZ <= alphas_moy)) {
//                 Lambda_min = Lambda_max;
//             } else {
//                 Lambda_max = Lambda_min;
//             }

//             if (fabs(Lambda_min - Lambda_max) <= 1e-5) {
//                 break;
//             }
//         }

//         return (Lambda_min + Lambda_max) / 2.0;
//     }

// double QCDParameters::computeAlphas(double Q, double Lambda, int nf) {
//     const double pi = 3.14159265358979323846;
//     double beta0 = 11.0 - 2.0 / 3.0 * nf;
//     double beta1 = 51.0 - 19.0 / 3.0 * nf;
//     double beta2 = 2857.0 - 5033.0 / 9.0 * nf + 325.0 / 27.0 * nf * nf;

//     double logRatio = log(pow(Q / Lambda, 2.0));
//     double logLogRatio = log(logRatio);
    
//     return 4.0 * pi / (beta0 * logRatio) * (1.0 - 2.0 * beta1 / (beta0 * beta0 * logRatio) * logLogRatio + 4.0 * beta1 * beta1 / (pow(beta0 * logRatio, 2.0)) * (pow(logLogRatio - 0.5, 2.0) + beta2 * beta0 / (8.0 * beta1 * beta1) - 5.0 / 4.0));
// }

double QCDParameters::runningAlphasCalculation(double Q, double Lambda, int nf) {
    if (Lambda5 == -1.0) {
        return -1.0;
    }
    double mb =4.19999981;
    double mtop = 172.399994;

    if (Lambda5 == 0.0) {
        Lambda = calculateLambda( mtop,  mb, 5);
        Lambda5 = Lambda;
    } else {
        Lambda = Lambda5;
    }

    if (Q <= mass_t && Q >= mass_b) {
        // 5 active flavors
        nf = 5;
    } else if (Q > mass_t) {
        // 6 active flavors
        nf = 6;
        Lambda = calculateLambda( mtop,  mb, nf);
    } else {
        // 4 active flavors
        nf = 4;
        Lambda = calculateLambda( mtop, mb, nf);
    }

    return calculateLambda(mtop, mb, nf);
}
