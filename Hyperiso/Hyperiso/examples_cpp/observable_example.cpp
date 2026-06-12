#include "Logger.h"
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "config.hpp"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    HyperisoMaster hyp; // Create the interface for hyperiso (reading/writing in the lha and all the options we want to use)
    HyperisoConfig config; // Config struct where we can put all the options we want for Hyperiso (general options)
    config.model = Model::SM; // The model we want to use, SM by default. If not THDM or SUSY, MARTY is needed.
    hyp.init("lha/si_input.flha", config); // Initialize program manager with LHA file and the config. Search in the Assets directory if relative path.
    LOG_INFO("HyperisoMaster initialized");

    QCDOrder order = QCDOrder::NNLO;

    ObservableInterface oi; // Initialize interface for observables calculation.

    oi.add_observable(Observables::BR_B__KSTAR_GAMMA, order); //An observable can be added with its name and the order at which the calculation will be performed.
    oi.add_observable(BinnedObservableId(Observables::F_L_B0__KSTAR0_MU_MU, {1.,6.}), order);//If the observable is binned, one can use the BinnedObservable class to precise the bin to be used for the calculation its also possible to use initializer list for simplicity.
    oi.add_observables(Decays::B__l_l, order); //If one want to calculate all observable of a certain decay, its possible to specify the decay using the add_observables API.

    //compute_observable return a vector of ObservableValue. An ObservableValue is composed of the id of the observable, the calculated value and the bin.
    LOG_INFO("F_L(B0 > K0* mu mu) (q²=[1, 6] GeV²) =", oi.compute_observable(Observables::F_L_B0__KSTAR0_MU_MU)[0].value);

    //One can use the get_current_observables API to get all observables available in the ObservableInterface.
    for (auto obs : oi.get_current_observables()) {
        LOG_INFO(obs);
    }

    //One can also get the experimental value and uncertainty of an observable (experimental values which are used in the Statistic part).
    LOG_INFO("Experimental value and uncertainty : of BR(B > K* gamma) =", oi.get_exp_value(Observables::BR_B__KSTAR_GAMMA), "+/-", oi.get_exp_uncertainty(Observables::BR_B__KSTAR_GAMMA, UncertaintyType::STAT));
    
    //If the value of an observable is present in multiple experiment (e.g CMS or LHCB) one can use the following API (using the struct ExperimentObs which contain a BinnedObservableId and the name of the experiment):
    LOG_INFO("Experimental value and uncertainty : of BR(B > K* gamma) =", oi.get_exp_value({"CMS",{Observables::F_L_B0__KSTAR0_MU_MU, {1., 6.}}}), "+/-", oi.get_exp_uncertainty({"CMS", {Observables::F_L_B0__KSTAR0_MU_MU, {1., 6.}}}, UncertaintyType::STAT));
    return 0;
}