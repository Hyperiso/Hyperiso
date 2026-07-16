#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "ParameterSetter.h"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

struct CsvRow {
    int step;
    double mu_b;
    std::string observable;
    std::optional<double> q2_min;
    std::optional<double> q2_max;
    double value;
};

int main() {
    HyperisoMaster hyp;
    HyperisoConfig config;
    config.model = Model::SM;
    hyp.init("lha/si_input.flha", config);

    QCDOrder order = QCDOrder::NNLO;
    ObservableInterface oi;
    ParameterSetter ps;

    oi.add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::A_T_2_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {1.1, 2}}, QCDOrder::NNLO, true)

      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::A_T_2_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {2, 4.3}}, QCDOrder::NNLO, true)

      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::A_T_2_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {4.3, 6}}, QCDOrder::NNLO, true)

      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::F_L_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::A_T_2_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_2_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_3_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_4_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_5_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_6_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true)
      .add_observable(BinnedObservableId{ObservableMapper::to_id(Observables::P_PRIME_8_B0__KSTAR0_MU_MU), {14.18, 16}}, QCDOrder::NNLO, true);

    std::vector<CsvRow> rows;
    rows.reserve(16 * 32);

    auto start = std::chrono::steady_clock::now();
    int step = 0;

    for (double mu_b = 1.0; mu_b < 5.0; mu_b += 0.05) {
        ++step;
        ps.mutate({ParameterType::WILSON, "B_SCALE", 1}, mu_b);
        const auto results = oi.compute_all();

        for (const auto& group : results) {
            const std::string obs_name = group.first.str();
            for (const auto& val : group.second) {
                CsvRow row;
                row.step = step;
                row.mu_b = mu_b;
                row.observable = obs_name;
                if (val.bin.has_value()) {
                    row.q2_min = val.bin->first;
                    row.q2_max = val.bin->second;
                }
                row.value = val.value;
                rows.push_back(std::move(row));
            }
        }
    }

    auto stop = std::chrono::steady_clock::now();
    const auto us = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();

    std::ofstream csv("hyperiso_bkstarll_scan.csv");
    csv << "step,mu_b,observable,q2_min,q2_max,value\n";
    csv << std::setprecision(17);
    for (const auto& row : rows) {
        csv << row.step << ','
            << row.mu_b << ','
            << row.observable << ',';
        if (row.q2_min.has_value()) {
            csv << *row.q2_min << ',' << *row.q2_max;
        } else {
            csv << ',';
        }
        csv << ',' << row.value << '\n';
    }
    csv.close();

    std::cout << "Temps calcul seul : " << us * 1e-6 << " s\n";
    std::cout << "num of steps : " << step << '\n';
    std::cout << "CSV : hyperiso_bkstarll_scan.csv\n";
    return 0;
}
