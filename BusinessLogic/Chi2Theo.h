#ifndef CHI2THEO_H
#define CHI2THEO_H

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <random>
#include <memory>
#include "Nuisance.h"
#include "ObservableStrategyChi2.h"
#include "Observable.h"
#include "../DataBase/json_parser.h"
#include "ObservableFactory.h"

class Chi2Theo {
public:
    static Chi2Theo& getInstance(const std::string& config_file = "");

    Chi2Theo(Chi2Theo const&) = delete;
    void operator=(Chi2Theo const&) = delete;

    void load_parameters_from_json(const std::string& config_file);

    // void calculate_observables();
    void print_observables() const;

    std::vector<std::vector<double>> calculate_covariance();

    std::map<std::string, Nuisance> parameters;
    std::vector<std::unique_ptr<Observable>> observables;
    std::vector<std::vector<double>> correlation_matrix;

private:
    Chi2Theo(const std::string& config_file);
    void initialize_observables();


};

#endif // CHI2THEO_H
