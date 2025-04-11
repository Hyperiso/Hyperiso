#include "MartyWilson.h"

void MartyWilson::LO_calculation() {

    std::unordered_set<ParamId> sources {
        {"EW_SCALE", 1}
    };

    auto func = [this] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src, std::shared_ptr<DependentParameter> dep_param) {
        double epsi = 1e-4;
        double ew_scale = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
        complex_t result = 0;
        for (size_t i = 0; i < df.getRowCount(); ++i) {
            double Q_match = df.iat<double>(i, "Q_match");
            if (fabs(Q_match - ew_scale) < epsi) {
                std::cout << this->get_name() << " waw" << std::endl;
                for (auto& _ : this->df.getColumnNames()) {
                    if (this->get_name()+"_real" == _) {
                        if (isnan(df.iat<double>(i, this->get_name()+"_real")) && isnan(df.iat<double>(i, this->get_name()+"_img"))) {
                            break;
                        }
                        std::cout << df.iat<double>(i, this->get_name()+"_real") << " BUTE" << std::endl;
                        result = {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")};
                        dep_param->set_expected(result);
                        return;
                        // this->set_WilsonCoeffMatching("LO", {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")});
                        //return {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")}; TODO
                    }
                } 
            }
        }
        MartyInterface MartyInterface;
        MartyInterface.calculate(this->get_name(), this->get_model(), /*TODO this->get_Q_match()*/ 81);
        df = csv_reader.read_csv(this->csv_path);
        df.setIndex(df.getColumn<double>("Q_match").to_string_vec());


        for (size_t i = 0; i < df.getRowCount(); ++i) {
            double Q_match = df.iat<double>(i, "Q_match");
            if (fabs(Q_match-/*TODO this->get_Q_match()*/ 0.) < epsi) {
                result = {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")};
                break;
                // this->set_WilsonCoeffMatching("LO", {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")});
                //return {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")}; TODO
            }
        }
        dep_param->set_expected(result);
    };

    WilsonParamComposer().compose_parameter(ParamId{"B_MATCH", LhaID(3040405, 4141, 2, 0)}, sources, func);

}