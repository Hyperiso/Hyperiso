#include "MartyWilson.h"


MartyWilson::MartyWilson(const LhaID& coeff_id, const std::string& storage_block, const std::string& model_name, const fs::path& model_path)
    : WilsonCoefficient(WCoefMapper::str(WCoefMapper::from_flha(coeff_id.get_parts()[0], coeff_id.get_parts()[1])), storage_block) {
    this->type = static_cast<ContributionType>(coeff_id.get_parts()[3]);
    this->set_model(model_name);
    std::unordered_set<ParamId> sources;

    std::string name = this->get_name();

    std::string csv_relative_path = "/MartyTemp/" + this->get_model() + "_wilson.csv";
    std::string csv_path = project_assets_root.data() +csv_relative_path;
    std::string marty_model = this->get_model();
    std::string marty_model_path = model_path;
    ContributionType cont = this->type;

    matching_info[QCDOrder::LO].compute = [&sources, name, csv_path, marty_model, marty_model_path] (const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) -> scalar_t {
        LOG_DEBUG("Updating coeff", name);
        double epsi = 1e-4;
        double ew_scale = src.at({ParameterType::WILSON, "EW_SCALE", 1})->get_val();
        scalar_t result;

        CSVReader csv_reader;
        DataFrame df;

        MartyInterface martyInterface;
        martyInterface.calculate(name, marty_model, ew_scale, marty_model_path);
        df = csv_reader.read_csv(csv_path);
        df.setIndex(df.getColumn<double>("Q_match").to_string_vec());

        for (size_t i = 0; i < df.getRowCount(); ++i) {
            double Q_match = df.iat<double>(i, "Q_match");
            if (fabs(Q_match - ew_scale) < epsi) {
                result = {df.iat<double>(i, name+"_real"), df.iat<double>(i, name+"_img")};
                break; 
            }
        }

        if (src.size() == 1) {
            std::set<std::string> special = martyInterface.get_special_blocks();
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
        }
        return result;
    };

    ParamId pid {ParameterType::WILSON, "EW_SCALE", 1};
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> dummy {{pid, std::make_shared<Parameter>(pid, 1, 0, 0)}};
    matching_info[QCDOrder::LO].compute(dummy);

    sources.emplace(ParamId{ParameterType::WILSON, "EW_SCALE", 1});
    matching_info[QCDOrder::LO].sources = sources;
    matching_info[QCDOrder::LO].lhaid = coeff_id;

    WCoef coef = WCoefMapper::from_flha(coeff_id.parts[0], coeff_id.parts[1]);
    ContributionType ct = static_cast<ContributionType>(coeff_id.parts[3]);
    matching_info[QCDOrder::NLO].lhaid = WCoefMapper::flha_full(coef, QCDOrder::NLO, ct);
    matching_info[QCDOrder::NNLO].lhaid = WCoefMapper::flha_full(coef, QCDOrder::NNLO, ct);
}