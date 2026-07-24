#include "Include.h"
#include "Logger.h"
#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "config.hpp"
#include "ParamOptimizerAdapter.h"


int main() {
        const std::string lha_path = "lha/si_input.flha";
        Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
        std::cerr << "[INFO] LHA = " << lha_path << "\n";

        HyperisoMaster hyp;
        HyperisoConfig config;
        config.model = Model::SM;

        hyp.init(lha_path, config);

        ObservableInterface oi;
        constexpr bool add_deps = false;

        oi.add_observables(Decays::M0_Mix, QCDOrder::LO);

        std::map<ObservableId, std::vector<ObservableValue>> all;
        all = oi.compute_all();

        std::ofstream out("mixing_obs.csv");
        out << "obs,value\n";
        out << std::setprecision(17);

        for (auto &[oid, value] : all) {
            out << ObservableMapper::str(oid) << "," << value[0].value << "\n";
        }
        
        out.close();

        return 0;
} 