#include <iostream>
#include "./Physical_Model/QCDParameters.h"

int main() {

    QCDParameters a;

    double lambda_1 = a.runningAlphasCalculation(1);
    double lambda_2 = a.runningAlphasCalculation(50);
    double lambda_3 = a.runningAlphasCalculation(100);
    double lambda_4 = a.runningAlphasCalculation(175);
    double lambda_5 = a.runningAlphasCalculation(201);
    double lambda_6 = a.runningAlphasCalculation(300);

    std::cout << lambda_1 << std::endl;
    std::cout << lambda_2 << std::endl;
    std::cout << lambda_3 << std::endl;
    std::cout << lambda_4 << std::endl;
    std::cout << lambda_5 << std::endl;
    std::cout << lambda_6 << std::endl;
    return 0;
}
