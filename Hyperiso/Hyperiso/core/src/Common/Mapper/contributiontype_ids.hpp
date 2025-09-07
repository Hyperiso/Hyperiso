
#pragma once
#include "generic_mapper.hpp"
#include "Map.h"

struct ContrTag {};
using ContributionTypeId = IdOf<ContrTag>;

class ContributionTypeMapper
: public GenericMapperNoExt<ContrTag, ContributionType, contributiontype_mapping>
{ public: using Base=GenericMapperNoExt<ContrTag,ContributionType,contributiontype_mapping>; using Base::str; };
