#include "Chi2Theo.h"
#include <fstream>

Chi2Theo::Chi2Theo(const std::string& config_file) {
    if (!config_file.empty()) {
        load_parameters_from_json(config_file);
    }
}

Chi2Theo& Chi2Theo::getInstance(const std::string& config_file) {
    static Chi2Theo instance(config_file);
    return instance;
}

void Chi2Theo::load_parameters_from_json(const std::string& config_file) {

    std::vector<Value> values;
    std::vector<Correlation> correlations;
    read_json(config_file, values, correlations);

    for (const auto& value : values) {
        parameters[value.name] = Nuisance(value.name, value.central_value, value.stat_error, value.syst_error);
    }

    for (const auto& obs : values) {
        std::vector<std::string> relevant_params = { obs.name };
        observables.push_back(ObservableFactory::createObservable("SpecificObservable", obs.name, relevant_params));
    }

    int n = values.size();
    correlation_matrix.resize(n, std::vector<double>(n, 0.0));

    std::unordered_map<std::string, int> param_index;
    int index = 0;
    for (const auto& value : values) {
        param_index[value.name] = index++;
    }

    for (const auto& correlation : correlations) {
        int idx1 = param_index[correlation.name1];
        int idx2 = param_index[correlation.name2];
        correlation_matrix[idx1][idx2] = correlation.value;
        correlation_matrix[idx2][idx1] = correlation.value;
    }
}

void Chi2Theo::calculate_observables() {

    for (auto& obs : observables) {
        obs->calculate(parameters);
    }
}

std::vector<std::vector<double>> Chi2Theo::calculate_covariance() {
    std::random_device rd;
    std::mt19937 rng(rd());
    int n = parameters.size();
    std::vector<std::vector<double>> covariance_matrix(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                covariance_matrix[i][j] += correlation_matrix[i][j] *
                                           (parameters.begin()->second.get_randomized_value(1, rng) - parameters.begin()->second.central_value) *
                                           (parameters.begin()->second.get_randomized_value(1, rng) - parameters.begin()->second.central_value);
            }
        }
    }
    return covariance_matrix;
}

void Chi2Theo::print_observables() const {
    for (const auto& obs : observables) {
        std::cout << "Observable: " << obs->name << ", Value: " << obs->value << std::endl;
    }
}