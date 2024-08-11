#include "chi2.h"
#include "Math.h"

void Chi2Manager::set_cov_tot() {

    for (auto& elem : chi_exp->get_covariance()) {
        cov_tot[elem.first] = elem.second + chi_theo->get_covariance()[elem.first];
        std::cout << elem.first.first << " is " << cov_tot[elem.first] <<  std::endl;
    }
}

void Chi2Manager::set_inv_cov_tot() {

    std::vector<std::string> diag_elem = getDiagonalElements(this->cov_tot);
    for (auto elem : diag_elem) {
        std::cout << "elem is " << elem << std::endl;
    }
    this->inv_cov_tot = invertMatrix(this->cov_tot, diag_elem);
}

void Chi2Manager::print_inv_cov() {

    for (auto& elem : this->inv_cov_tot) {
        std::cout << "[" << elem.first.first << " , " << elem.first.second << "]" << " : " << elem.second << std::endl;
    }
}