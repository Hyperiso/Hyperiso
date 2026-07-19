#include "MartyWilson.h"


MartyWilson::MartyWilson(MartyWilsonConfig config)
    : WilsonCoefficient(WCoefMapper::str(WCoefMapper::from_flha(config.coeff_id.get_parts()[0], config.coeff_id.get_parts()[1])), config.storage_block) {
    if (config.model_path.empty()) {
        throw std::invalid_argument(
            "MartyWilson requires an explicit runtime-resolved MARTY model path"
        );
    }
    this->type = static_cast<ContributionType>(config.coeff_id.get_parts()[3]);
    this->set_model(config.model_name);
    std::string name = this->get_name();

    std::string csv_path = config.csv_path;
    std::string marty_model = config.model_name;
    std::string marty_generation_model = config.generation_model_name;
    bool marty_sm_like_filter = config.sm_like_filter;
    bool marty_bsm_split_generation = config.bsm_split_generation;
    bool marty_full_target_generation = config.full_target_generation;
    std::string marty_model_path = config.model_path;

    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy = config.marty_proxy;
    ContributionType contribution_type = this->type;

    matching_info[QCDOrder::LO].compute = [name, csv_path, marty_model, marty_generation_model, marty_sm_like_filter, marty_bsm_split_generation, marty_full_target_generation, marty_model_path, marty_proxy, contribution_type] (const ParamSrc& src) -> scalar_t {
        LOG_DEBUG("Updating coeff", name);
        double epsi = 1e-4;
        double ew_scale = src.get_val({ParameterType::WILSON, "EW_SCALE", 1});
        if (!marty_proxy) {
            throw std::runtime_error(
                "MartyWilson cannot evaluate '" + name + "': no MARTY proxy was configured"
            );
        }

        scalar_t result{};
        bool found_matching_scale = false;

        CSVReader csv_reader;
        DataFrame df;

        marty_proxy->calculate(
            name,
            marty_model,
            marty_generation_model,
            ew_scale,
            marty_model_path,
            marty_sm_like_filter,
            marty_bsm_split_generation,
            marty_full_target_generation
        );
        df = csv_reader.read_csv(csv_path);
        df.setIndex(df.getColumn<double>("Q_match").to_string_vec());

        std::string csv_column_base = name;
        if (marty_bsm_split_generation && name == "CP10") {
            if (contribution_type == ContributionType::SM) {
                csv_column_base = name + "_SM_SPLIT";
            } else if (contribution_type == ContributionType::BSM) {
                csv_column_base = name + "_BSM_SPLIT";
            } else if (contribution_type == ContributionType::TOTAL) {
                csv_column_base = name + "_TOTAL_SPLIT";
            }
        }

        for (size_t i = 0; i < df.getRowCount(); ++i) {
            double Q_match = df.iat<double>(i, "Q_match");
            if (fabs(Q_match - ew_scale) < epsi) {
                result = {df.iat<double>(i, csv_column_base+"_real"), df.iat<double>(i, csv_column_base+"_img")};
                found_matching_scale = true;
                break;
            }
        }

        if (!found_matching_scale) {
            throw std::runtime_error(
                "MartyWilson CSV for '" + name + "' contains no row at Q_match="
                + std::to_string(ew_scale)
            );
        }

        return result;
    };

    ParamId pid {ParameterType::WILSON, "EW_SCALE", 1};
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> dummy {{pid, std::make_shared<Parameter>(pid, 1, 0, 0)}};
    matching_info[QCDOrder::LO].compute(ParamSrc(dummy));

    std::unordered_set<ParamId> sources;
    const std::set<std::string> special = marty_proxy->get_special_blocks();
    for (const auto& par : marty_proxy->get_dependencies(name)) {
        if (special.contains(par.block)) {
            continue;
        }
        sources.emplace(ParamId{
            par.is_bsm ? ParameterType::BSM : ParameterType::SM,
            par.block,
            par.code
        });
    }
    sources.emplace(ParamId{ParameterType::WILSON, "EW_SCALE", 1});
    matching_info[QCDOrder::LO].sources = std::move(sources);
    matching_info[QCDOrder::LO].lhaid = config.coeff_id;

    WCoef coef = WCoefMapper::from_flha(config.coeff_id.parts[0], config.coeff_id.parts[1]);
    ContributionType ct = static_cast<ContributionType>(config.coeff_id.parts[3]);
    matching_info[QCDOrder::NLO].lhaid = WCoefMapper::flha_full(coef, QCDOrder::NLO, ct);
    matching_info[QCDOrder::NNLO].lhaid = WCoefMapper::flha_full(coef, QCDOrder::NNLO, ct);
}