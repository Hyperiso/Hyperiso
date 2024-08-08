#ifndef CHI2EXP_H
#define CHI2EXP_H

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


class Chi2Exp {
public:
    static Chi2Exp& getInstance(const std::string& config_file = "");

    Chi2Exp(Chi2Exp const&) = delete;
    void operator=(Chi2Exp const&) = delete;


    // void calculate_observables();
    void print_observables() const;
    void print_correlations() const;
    void print_correlations_matrix();
    void fill_from_theory();

    void calculate_covariance();

    // std::map<std::string, Nuisance> parameters;
    std::vector<std::unique_ptr<Observable>> observables;
    std::map<std::pair<std::string, std::string>, double> correlation_matrix;

private:
    Chi2Exp(const std::string& config_file);
    void initialize_observables(const std::string& config_file);

    std::vector<Value> values;
    std::vector<Correlation> correlations;

    ObservableMapper& mapper = ObservableMapper::getInstance();
};

#endif // CHI2EXP_H
