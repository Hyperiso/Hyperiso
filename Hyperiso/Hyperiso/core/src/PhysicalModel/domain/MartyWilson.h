#ifndef MARTY_WILSON_H
#define MARTY_WILSON_H

#include <complex>
#include "Wilson.h"
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

    MartyWilson(MartyWilsonConfig config);

    std::string get_model() {
        return this->model;
    }
    void set_model(std::string model) {
        this->model = model;
    }

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<MartyWilson>(*this);
    }

private:
    std::string model{"SM"};

};

#endif