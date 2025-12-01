#ifndef MAP_H
#define MAP_H

#include <map>
#include <string>
#include <vector>

#include "GeneralEnum.h"
#include "General.h"

const std::map<Observables, std::string>& observable_mapping();
const std::map<Observables, LhaID>&       observable_flha_mapping();

const std::map<QCDOrder,        std::string>& order_mapping();
const std::map<WGroup,          std::string>& group_mapping();
const std::map<WCoef,           std::string>& wcoef_mapping();
const std::map<WCoef, std::pair<int,int>>&    wcoef_flha_mapping();
const std::map<ParameterType,   std::string>& parametertype_mapping();
const std::map<Model,           std::string>& model_mapping();
const std::map<WilsonBasis,     std::string>& wilsonbasis_mapping();
const std::map<ContributionType,std::string>& contributiontype_mapping();
const std::map<MassType,        std::string>& masstype_mapping();
const std::map<ScaleType,       std::string>& scaletype_mapping();
const std::map<Decays,          std::string>& decays_mapping();
const std::map<Decays, std::vector<Observables>>& decay_observable_mapping();

#endif // MAP_H
