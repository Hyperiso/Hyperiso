#include "Chi2Theo.h"
#include <fstream>

Chi2Theo::Chi2Theo(const std::string& config_file) {
    std::cout << "Creating Chi2Theo" << std::endl;
    if (!config_file.empty()) {
        load_parameters_from_json(config_file);
        std::cout << "JSON loaded" << std::endl;
    }
    std::cout << "trying to initialize obs" << std::endl;
    initialize_observables();
    std::cout << "observable initialized" << std::endl;

}

Chi2Theo* Chi2Theo::GetInstance(const std::string& config_file) {
        if (!Chi2Theo::instance) {
            Chi2Theo::instance = new Chi2Theo(config_file);
        }
        return Chi2Theo::instance;
    }

void Chi2Theo::load_parameters_from_json(const std::string& config_file) {
    std::vector<Value> values;
    std::vector<Correlation> correlations;
    read_json(config_file, values, correlations);

    for (const auto& value : values) {
        parameters[value.name] = *ObservableFactory::createNuisance(value.name, value.central_value, value.stat_error, value.syst_error);
    }

    int n = values.size();

    std::unordered_map<std::string, int> param_index;
    int index = 0;
    for (const auto& value : values) {
        param_index[value.name] = index++;
    }

    for (const auto& correlation : correlations) {
        int idx1 = param_index[correlation.name1];
        int idx2 = param_index[correlation.name2];
        // correlation_matrix[idx1][idx2] = correlation.value;
        // correlation_matrix[idx2][idx1] = correlation.value;
    }
}

// void Chi2Theo::calculate_observables() {

//     for (auto& obs : observables) {
//         obs->evaluate();
//     }
// }

std::map<std::pair<std::string, std::string>, double> Chi2Theo::calculate_covariance() {
    std::random_device rd;
    std::mt19937 rng(rd());
    int n = parameters.size();
    std::map<std::pair<std::string, std::string>, double> covariance_matrix;

    // for (int i = 0; i < n; ++i) {
    //     for (int j = 0; j < n; ++j) {
    //         for (int k = 0; k < n; ++k) {
    //             covariance_matrix[i][j] += correlation_matrix[i][j] *
    //                                        (parameters.begin()->second.get_randomized_value(1, rng) - parameters.begin()->second.central_value) *
    //                                        (parameters.begin()->second.get_randomized_value(1, rng) - parameters.begin()->second.central_value);
    //         }
    //     }
    // }
    return covariance_matrix;
}

void Chi2Theo::print_observables() const {
    std::cout << "trying to print observables" << std::endl;
    for (const auto& obs : observables) {
        std::cout << "Observable: " << static_cast<int>(obs->getId()) << ", Value: " << obs->getValue() << std::endl;
    }
}

void Chi2Theo::initialize_observables() {
    for (int i = static_cast<int>(Observables::FIRST)+1; i < static_cast<int>(Observables::LAST); ++i) {
        std::cout << "new observable" << std::endl;
        Observables obs_id = static_cast<Observables>(i);
        auto observable = ObservableFactory::createObservable(obs_id, 0, 2, 5.27958, 0);
        observables.push_back(std::move(observable));
        std::cout << "new observable ended" << std::endl;
    }
    std::cout << "end observable" << std::endl;
}

Chi2Theo* Chi2Theo::instance = nullptr;