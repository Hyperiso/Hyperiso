#include "Chi2Exp.h"

Chi2Exp::Chi2Exp(const std::string& config_file) {
    std::cout << "Creating Chi2Theo" << std::endl;
    if (!config_file.empty()) {
        std::cout << "trying to initialize obs" << std::endl;
        initialize_observables(config_file);
        std::cout << "observable initialized" << std::endl;
        std::cout << "JSON loaded" << std::endl;
    }
    

}

Chi2Exp& Chi2Exp::getInstance(const std::string& config_file) {
    static Chi2Exp instance(config_file);
    return instance;
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
    }
    
    correlation_matrix = temp_step;
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

    calculate_covariance();
    for (const auto& truc : correlation_matrix) {
        std::cout << "With pair format, " <<"between " <<  truc.first.first << " and " << truc.first.second << " : " << truc.second << std::endl;
    }

}