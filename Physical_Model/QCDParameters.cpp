#include "QCDParameters.h"

double QCDParameters::calculateLambda(double MZ, double alphas_MZ, int nf) {
        const double pi = 3.14159265358979323846;
        double beta0 = 11.0 - (2.0 / 3.0) * nf;
        double beta1 = 51.0 - (19.0 / 3.0) * nf;
        double beta2 = 2857.0 - (5033.0 / 9.0) * nf + (325.0 / 27.0) * nf * nf;

        double Lambda_min = 1e-3;
        double Lambda_max = 1.0;
        double alphas_min, alphas_moy;

        while (true) {
            double log_MZ_Lambda_min = log(pow(MZ / Lambda_min, 2));
            alphas_min = 4.0 * pi / (beta0 * log_MZ_Lambda_min) *
                         (1.0 - 2.0 * beta1 / (beta0 * beta0) * log(log_MZ_Lambda_min) / log_MZ_Lambda_min +
                          4.0 * beta1 * beta1 / (beta0 * beta0 * beta0 * beta0 * log_MZ_Lambda_min * log_MZ_Lambda_min) *
                          (pow(log(log_MZ_Lambda_min) - 0.5, 2) + beta2 * beta0 / (8.0 * beta1 * beta1) - 5.0 / 4.0));

            double Lambda_moy = (Lambda_min + Lambda_max) / 2.0;
            double log_MZ_Lambda_moy = log(pow(MZ / Lambda_moy, 2));
            alphas_moy = 4.0 * pi / (beta0 * log_MZ_Lambda_moy) *
                         (1.0 - 2.0 beta1 / (beta0 * beta0) * log(log_MZ_Lambda_moy) / log_MZ_Lambda_moy +
                          4.0 * beta1 * beta1 / (beta0 * beta0 * beta0 * beta0 * log_MZ_Lambda_moy * log_MZ_Lambda_moy) *
                          (pow(log(log_MZ_Lambda_moy) - 0.5, 2) + beta2 * beta0 / (8.0 * beta1 * beta1) - 5.0 / 4.0));

            if ((alphas_MZ >= alphas_min) && (alphas_MZ <= alphas_moy)) {
                Lambda_max = Lambda_moy;
            } else {
                Lambda_min = Lambda_moy;
            }

            if (fabs(Lambda_min - Lambda_max) <= 1e-5) {
                break;
            }
        }

        return (Lambda_min + Lambda_max) / 2.0;
    }


double QCDParameters::calculateLambda(double MZ, double alphas_MZ, int nf) {
        const double pi = 3.14159265358979323846;
        double beta0 = 11.0 - (2.0 / 3.0) * nf;
        double beta1 = 51.0 - (19.0 / 3.0) * nf;
        double beta2 = 2857.0 - (5033.0 / 9.0) * nf + (325.0 / 27.0) * nf * nf;

        double Lambda_min = 1e-3;
        double Lambda_max = 1.0;
        double alphas_min, alphas_moy;

        while (true) {
            double log_MZ_Lambda_min = log(pow(MZ / Lambda_min, 2));
            alphas_min = 4.0 * pi / (beta0 * log_MZ_Lambda_min) *
                         (1.0 - 2.0 * beta1 / (beta0 * beta0) * log(log_MZ_Lambda_min) / log_MZ_Lambda_min +
                          4.0 * beta1 * beta1 / (beta0 * beta0 * beta0 * beta0 * log_MZ_Lambda_min * log_MZ_Lambda_min) *
                          (pow(log(log_MZ_Lambda_min) - 0.5, 2) + beta2 * beta0 / (8.0 * beta1 * beta1) - 5.0 / 4.0));

            Lambda_max = (Lambda_min + Lambda_max) / 2.0;
            double log_MZ_Lambda_max = log(pow(MZ / Lambda_max, 2));
            alphas_moy = 4.0 * pi / (beta0 * log_MZ_Lambda_max) *
                         (1.0 - 2.0 * beta1 / (beta0 * beta0) * log(log_MZ_Lambda_max) / log_MZ_Lambda_max +
                          4.0 * beta1 * beta1 / (beta0 * beta0 * beta0 * beta0 * log_MZ_Lambda_max * log_MZ_Lambda_max) *
                          (pow(log(log_MZ_Lambda_max) - 0.5, 2) + beta2 * beta0 / (8.0 * beta1 * beta1) - 5.0 / 4.0));

            if ((alphas_MZ >= alphas_min) && (alphas_MZ <= alphas_moy)) {
                Lambda_min = Lambda_max;
            } else {
                Lambda_max = Lambda_min;
            }

            if (fabs(Lambda_min - Lambda_max) <= 1e-5) {
                break;
            }
        }

        return (Lambda_min + Lambda_max) / 2.0;
    }

double QCDParameters::runningAlphasCalculation(double Q, double Lambda, int nf) {
        const double pi = 3.14159265358979323846;
        double beta0 = 11.0 - (2.0 / 3.0) * nf;
        double beta1 = 51.0 - (19.0 / 3.0) * nf;
        double beta2 = 2857.0 - (5033.0 / 9.0) * nf + (325.0 / 27.0) * nf * nf;

        double logQ2Lambda2 = log(pow(Q / Lambda, 2));
        double alphas = 4.0 * pi / (beta0 * logQ2Lambda2) *
                        (1.0 - 2.0 * beta1 / (beta0 * beta0) * log(logQ2Lambda2) / logQ2Lambda2 +
                         4.0 * beta1 * beta1 / (beta0 * beta0 * beta0 * beta0 * logQ2Lambda2 * logQ2Lambda2) *
                         (pow(log(logQ2Lambda2) - 0.5, 2) + beta2 * beta0 / (8.0 * beta1 * beta1) - 5.0 / 4.0));

        return alphas;
    }
