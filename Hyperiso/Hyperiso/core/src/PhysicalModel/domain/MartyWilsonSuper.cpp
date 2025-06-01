#include "MartyWilsonSuper.h"


matching_info[QCDOrder::NNLO] = {
        {
            {"WPARAM_MATCH_SM", 3},
            {"WPARAM_MATCH_SM", LhaID(2, 1)}
        },
        [](const auto& src) {
            auto L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
            auto xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
            return -T(xt) + 7987./72. + 17. * PI2 / 3. + 475./6. * L + 17. * L * L;
        },
        LhaID(3040405, 6161, 2, 0)
    };


MartyWilson::MartyWilson(const std::string& coeff_name, ContributionType cont_type)
    : WilsonCoefficient(coeff_name, "B_MATCH") {
    this->set_name(coeff_name);
    this->type = cont_type;

    matching_info[QCDOrder::LO] = {
        {
            {"WPARAM_MATCH_SM", 3},
            {"WPARAM_MATCH_SM", LhaID(2, 1)}
        },
        [](const auto& src) {
            auto L = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 3})->get_val();
            auto xt = src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {2, 1}})->get_val();
            return -T(xt) + 7987./72. + 17. * PI2 / 3. + 475./6. * L + 17. * L * L;
        },
        LhaID(3040405, 6161, 2, 0)
    };

}

void MartyWilson::LO_calculation() {

    LOG_INFO("In MartyWilson::LO_calculation of coefficient " + this->get_name());

    std::unordered_set<ParamId> sources {
        {"EW_SCALE", 1}
    };

    std::string name = this->get_name();

    std::string csv_relative_path = "/MartyTemp/" + this->get_model() + "_wilson.csv";
    std::string csv_path = project_assets_root.data() +csv_relative_path;
    std::string mod = this->get_model();
    ContributionType cont = this->type;
    auto func = [&sources, name, csv_path, mod, cont] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        LOG_INFO("Updating coeff");
        double epsi = 1e-4;
        double ew_scale = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
        complex_t result = 0;

        CSVReader csv_reader;
        DataFrame df;

        // for (size_t i = 0; i < df.getRowCount(); ++i) {
        //     double Q_match = df.iat<double>(i, "Q_match");
        //     if (fabs(Q_match - ew_scale) < epsi) {
        //         std::cout << this->get_name() << " waw" << std::endl;
        //         for (auto& _ : this->df.getColumnNames()) {
        //             if (this->get_name()+"_real" == _) {
        //                 if (isnan(df.iat<double>(i, this->get_name()+"_real")) && isnan(df.iat<double>(i, this->get_name()+"_img"))) {
        //                     break;
        //                 }
        //                 std::cout << df.iat<double>(i, this->get_name()+"_real") << " BUTE" << std::endl;
        //                 result = {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")};
        //                 dep_param->set_expected(result);
        //                 return;
        //                 // this->set_WilsonCoeffMatching("LO", {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")});
        //                 //return {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")}; TODO
        //             }
        //         } 
        //     }
        // }

        MartyInterface martyInterface;
        martyInterface.calculate(name, mod, ew_scale);
        df = csv_reader.read_csv(csv_path);
        df.setIndex(df.getColumn<double>("Q_match").to_string_vec());

        for (size_t i = 0; i < df.getRowCount(); ++i) {
            double Q_match = df.iat<double>(i, "Q_match");
            if (fabs(Q_match - ew_scale) < epsi) {
                result = {df.iat<double>(i, name+"_real"), df.iat<double>(i, name+"_img")};
                break;
                // this->set_WilsonCoeffMatching("LO", {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")});
                //return {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")}; 
            }
        }

        std::set<std::string> special = {"KIN", "WEIN", "Finite", "REGPROP"}; //TODO : do better with this, in SMParamSetter
        for (auto &par : martyInterface.get_dependencies(name)) {
            if (std::find(special.begin(), special.end(),par.block) != special.end()) {
                continue;
            }
            if (par.is_bsm) {
                sources.emplace(ParamId{ParameterType::BSM, par.block, par.code});
            } else {
                sources.emplace(ParamId{ParameterType::SM, par.block, par.code});
            }
        }

        dep_param->set_expected(result);
    };
    ParamId pid {ParameterType::WILSON, "EW_SCALE", 1};
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> dummy {{pid, std::make_shared<Parameter>(pid, 1, 0, 0)}};
    func(dummy, std::make_shared<DependentParameter>(dummy, func));

    //TODO  :configuration here to SM, need to put at FULL but cannot work with WilsonGroup at the moment (only sm)
    WilsonParamComposer().compose_parameter(ParamId{this->storage_block, WCoefMapper::flha_full(WCoefMapper::enum_elt(this->coeffName), QCDOrder::LO, this->type)}, sources, func);

}