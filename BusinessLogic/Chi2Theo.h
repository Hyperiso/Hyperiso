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
    static Chi2Theo* instance;
public:
    static Chi2Theo* GetInstance(const std::string& config_file = "");

    Chi2Theo(Chi2Theo const&) = delete;
    void operator=(Chi2Theo const&) = delete;

    void load_parameters_from_json(const std::string& config_file);

    std::map<std::string, std::complex<double>> get_obs() {return obs;}
    std::map<std::pair<std::string, std::string>, double> get_covariance() {return correlation_matrix;}
    // void calculate_observables();
    void print_observables() const;

    std::map<std::string, Nuisance> parameters;
    void calculate_covariance();

    std::vector<std::unique_ptr<Observable>> observables;

    std::map<std::string, std::complex<double>> obs;
    std::map<std::pair<std::string, std::string>, double> correlation_matrix;

private:
    Chi2Theo(const std::string& config_file);
    void initialize_observables();

    ObservableMapper* mapper = ObservableMapper::GetInstance();
    void transform_obs();
};

#endif // CHI2THEO_H
