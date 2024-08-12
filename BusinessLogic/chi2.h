#include "Chi2Exp.h"
#include "Chi2Theo.h"

class Chi2Manager {

    int model;
    int order;
    double scale;
    int wilson_basis;

    Chi2Exp* chi_exp = Chi2Exp::GetInstance("../../DataBase/data_exp.json");
    Chi2Theo* chi_theo = Chi2Theo::GetInstance("../../DataBase/data_theo.json");
public:
    Chi2Manager(int model, int order, double scale, int wilson_basis) :
    model(model), order(order), scale(scale), wilson_basis(wilson_basis) {
        set_cov_tot();
        set_inv_cov_tot();
    }

    std::map<std::pair<std::string, std::string>, double> cov_tot;
    std::map<std::pair<std::string, std::string>, double> inv_cov_tot;

    void set_cov_tot();
    void set_inv_cov_tot();

    double get_chi2();

    void print_inv_cov();
};