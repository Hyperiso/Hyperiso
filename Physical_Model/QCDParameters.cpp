#include "./QCDParameters.h"
#include <iostream>

QCDParameters::QCDParameters() {
    mass_Z = 91.1699982;
    alphas_MZ = 0.117200002;
    
}

double QCDParameters::alphasRunning(double Q, double Lambda, int nf) const{
    double beta0 = 11.-2./3.*nf;
    double beta1=51.-19./3.*nf;
    double beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

    return 4.*pi/beta0/log(pow(Q/Lambda,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda,2.)))/log(pow(Q/Lambda,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda,2.)),2.)*(pow(log(log(pow(Q/Lambda,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
}

double QCDParameters::DichotomieLambda(double alpha_running, double Q, int nf){

    double Lambda_min=1.e-3;
    double Lambda_max=1.;
    double alphas_min=0.;
    double alphas_max=0;
    double alphas_moy = 0;
    double Lambda_moy = 0;
    while((fabs(1.-alphas_min/alpha_running)>=1.e-4)&&(fabs(1.-Lambda_min/Lambda_max)>1.e-5))
    {
        alphas_min = alphasRunning(Q,Lambda_min, nf);
        alphas_max = alphasRunning(Q,Lambda_max, nf);

        Lambda_moy=(Lambda_min+Lambda_max)/2.;
        alphas_moy = alphasRunning(Q,Lambda_moy, nf);

        if((alpha_running>=alphas_min)&&(alpha_running<=alphas_moy))
            Lambda_max=Lambda_moy;
        else Lambda_min=Lambda_moy;
    }
    
    if(fabs(1.-Lambda_min/Lambda_max)<=1.e-5)
    {
        Lambda5=-1.;
        return -1.;
    }
    return Lambda_min;
}

double QCDParameters::runningAlphasCalculation(double Q){
    if (Lambda5 == -1.0) {
        return -1.0;
    }

    if (Lambda5 == 0.0) {
        Lambda5 = DichotomieLambda( alphas_MZ,  mass_Z, 5);
        if (Lambda5 == -1);
            return Lambda5;
    }
    double alphas_running = alphasRunning(Q, Lambda5, 5);

    if (Q <= mass_t && Q >= mass_b) {
        // 5 active flavors
        return alphas_running;

    } else if (Q > mass_t) {
        // 6 active flavors
        nf = 6;
        alphas_running = alphasRunning(mass_t, Lambda5, 5);
        Lambda6 = DichotomieLambda( alphas_running,  mass_t, nf);
        return alphasRunning(Q,Lambda6, nf);
    } else {
        // 4 active flavors
        alphas_running = alphasRunning(mass_b, Lambda5, 5);
        nf = 4;
        Lambda4 = DichotomieLambda( alphas_running,  mass_b, nf);
        return alphasRunning(Q,Lambda4,  nf);
    }

    throw std::logic_error("Invalid parameters or conditions in alphas_running function");
}
