#include <iostream>
#include "./Physical_Model/QCDParameters.h"

int main() {

    QCDParameters a;

    double mb =4.4;
    double mtop = 170;

    double lambda_1 = a.calculateLambda(mtop, mb,1);
    double lambda_2 = a.calculateLambda(mtop, mb,2);
    double lambda_3 = a.calculateLambda(mtop, mb,3);
    double lambda_4 = a.calculateLambda(mtop, mb,4);
    double lambda_5 = a.calculateLambda(mtop, mb,5);
    double lambda_6 = a.calculateLambda(mtop, mb,6);

    std::cout << lambda_1 << std::endl;
    std::cout << lambda_2 << std::endl;
    std::cout << lambda_3 << std::endl;
    std::cout << lambda_4 << std::endl;
    std::cout << lambda_5 << std::endl;
    std::cout << lambda_6 << std::endl;
    return 0;
}
