#pragma once
#include <complex>
#include "Wilson.h"
// #include "Interpolator.h"
#include "config.hpp"
#include "DataFrame.h"
#include "CSVReader.h"
#include "InterpretedParam.h" //Only for template argument !!
#include "IMartyWilsonProxy.h"
#include "config.hpp"
#include "Utils.h"
#include <iostream>
#include <math.h>


struct MartyWilsonConfig {
    std::string model_name{"SM"};
    fs::path model_path{project_assets_root.data() + std::string() + "/input_files/marty_model/sm.h"};
    std::string csv_path{project_assets_root.data() + std::string() + +"/MartyTemp/SM_wilson.csv"};
    LhaID coeff_id;
    std::string storage_block;
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy;

    MartyWilsonConfig(const LhaID& id, const std::string& storage_block_name, fs::path model_path, std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> proxy)
        : coeff_id(id), storage_block(storage_block_name), model_path(model_path), marty_proxy(proxy) {}

    MartyWilsonConfig(const std::string& model_name,const LhaID& id, const std::string& storage_block_name, fs::path model_path, std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> proxy)
        : model_name(model_name), coeff_id(id), storage_block(storage_block_name), model_path(model_path), marty_proxy(proxy) {
            csv_path = project_assets_root.data() + std::string() +"/MartyTemp/" + model_name + "_wilson.csv";
        }
};


class MartyWilson : public WilsonCoefficient {
public:
    // MartyWilson(const std::string& coeff_name, const std::string& storage_block)
    //     : WilsonCoefficient(coeff_name, storage_block) {
    //     this->type = ContributionType::TOTAL;
    // }

    MartyWilson(MartyWilsonConfig config);

    std::string get_model() {
        return this->model;
    }
    void set_model(std::string model) {
        this->model = model;
    }

    // void LO_calculation();

    // void NLO_calculation() override {} //TODO, at least deal properly
    // void NNLO_calculation() override {} //TODO

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<MartyWilson>(*this);
    }

private:
    std::string model{"SM"};

    // std::pair<size_t, size_t> find_closest_Q_matches(double target_Q_match) {
        
    //     size_t closest_below = 0, closest_above = 0;
    //     double min_diff_below = std::numeric_limits<double>::max();
    //     double min_diff_above = std::numeric_limits<double>::max();

    //     for (size_t i = 0; i < df.getRowCount(); ++i) {
    //         double Q_match = df.iat<double>(i, "Q_match");
    //         double diff = Q_match - target_Q_match;
    //         std::cout << "Q_match : " << Q_match << std::endl;
    //         if (diff <= 0 && std::abs(diff) < min_diff_below) {
    //             closest_below = i;
    //             min_diff_below = std::abs(diff);
    //         } else if (diff > 0 && diff < min_diff_above) {
    //             closest_above = i;
    //             min_diff_above = diff;
    //         }
    //     }
    //     std::cout << closest_above << " " << closest_below << std::endl;
    //     return {closest_below, closest_above};
    // }
};