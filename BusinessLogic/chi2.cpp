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

double Chi2Manager::get_chi2() {

    double result{0};
    for (auto& obs1 : chi_theo->get_obs()) {
        for (auto& obs2 : chi_theo->get_obs()) {
            std::cout << "real of theory part : " << std::real(obs1.second) << std::endl;
            result += (std::real(obs1.second)-std::real(chi_exp->get_obs()[obs1.first])) * this->inv_cov_tot[std::make_pair(obs1.first, obs2.first)] * (std::real(obs2.second)-std::real(this->chi_exp->get_obs()[obs2.first]));

        }
    }
    return result;
}