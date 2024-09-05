#include "Chi2Exp.h"
ObservableMapper* ObservableMapper::instance = nullptr;

Chi2Exp::Chi2Exp(const std::string& config_file) {
    std::cout << "Creating Chi2Exp" << std::endl;
    if (!config_file.empty()) {
        std::cout << "trying to initialize obs" << std::endl;
        initialize_observables(config_file);
        std::cout << "observable initialized" << std::endl;
        std::cout << "JSON loaded" << std::endl;
    }
    calculate_covariance();
    fill_from_theory();
}

Chi2Exp* Chi2Exp::GetInstance(const std::string& config_file) {
        if (!Chi2Exp::instance) {
            Chi2Exp::instance = new Chi2Exp(config_file);
        }
        return Chi2Exp::instance;
    }

void Chi2Exp::initialize_observables(const std::string& config_file) {

    
    read_json(config_file, values, correlations);
    
}

void Chi2Exp::calculate_covariance() {

    std::map<std::pair<std::string, std::string>, double> temp_step{};

    for (const auto& value : values) {
        temp_step[std::make_pair(value.name, value.name)] = 1;
    }

    for (const auto& correlation : correlations) {
        temp_step[std::make_pair(correlation.name1,correlation.name2)] = correlation.value;
        temp_step[std::make_pair(correlation.name2,correlation.name1)] = correlation.value;
    }
    
    correlation_matrix = temp_step;
}

void Chi2Exp::fill_from_theory() {

    for (const auto& obs : mapper->get_map()) {
        if (!correlation_matrix.contains(std::make_pair(obs.second, obs.second)) ){
            std::cout << "yes yes " << obs.second << std::endl;
            correlation_matrix[std::make_pair(obs.second, obs.second)] = 1;
        }
        else {
            std::cout << "no noo " << obs.second << std::endl;
        }
    }
}
void Chi2Exp::print_observables() const{

    for (const auto& value : values) {
        std::cout << value.name << " : " << value.central_value << ", " << value.stat_error << " + " << value.syst_error << std::endl;
    }
}

void Chi2Exp::print_correlations() const{

    for (const auto& correlation : correlations) {
        std::cout << "Between " <<  correlation.name1 << " and " << correlation.name2 << " : " << correlation.value << std::endl;
    }
}

void Chi2Exp::print_correlations_matrix() {

    
    for (const auto& truc : correlation_matrix) {
        std::cout << "With pair format, " <<"between " <<  truc.first.first << " and " << truc.first.second << " : " << truc.second << std::endl;
    }

}

Chi2Exp* Chi2Exp::instance = nullptr;