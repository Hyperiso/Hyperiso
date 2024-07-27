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

void Chi2Exp::print_observables() const{

    for (const auto& value : values) {
        std::cout << value.name << " : " << value.central_value << ", " << value.stat_error << " + " << value.syst_error << std::endl;
    }
}