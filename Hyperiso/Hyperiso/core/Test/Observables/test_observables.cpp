#include "Logger.h"
#include "ObservablesUtils.h"
#include "config.hpp"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    const std::string lha_path = "lha/si_input.flha";
    const Model model = Model::SM; 
    const QCDOrder order = QCDOrder::NNLO; 
    const std::string out_dir = project_root.data()+ std::string("Test/csv/observables");
    const std::string ref_dir = project_root.data()+ std::string("Test/csv/superiso/observables");
    const double tolerance = 1e-4;

    runObservablesTest(lha_path, model, order, out_dir, ref_dir, tolerance);
    return 0;
}
