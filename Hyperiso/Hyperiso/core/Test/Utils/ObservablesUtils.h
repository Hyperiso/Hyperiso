#ifndef OBSERVABLES_UTILS_H
#define OBSERVABLES_UTILS_H

#include <string>
#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "config.hpp"

void writeDecayObservablesCsv(ObservableInterface& oi,
                              Decays dec,
                              QCDOrder order,
                              const std::string& out_dir);

void writeAllDecaysObservablesCsv(ObservableInterface& oi,
                                  QCDOrder order,
                                  const std::string& out_dir);

void runObservablesTest(const std::string& lha_path,
                        Model model,
                        QCDOrder order,
                        const std::string& out_dir,
                        const std::string& ref_dir,
                        double tolerance);

#endif
