#ifndef CONTRIBUTION_TYPE_IDS_H
#define CONTRIBUTION_TYPE_IDS_H

#include "generic_mapper.h"
#include "Map.h"

/**
 * @file ContributionTypeMapper.h
 * @brief Mapper for ContributionType <-> string identifiers.
 */
struct ContrTag {};
using ContributionTypeId = IdOf<ContrTag>;

class ContributionTypeMapper
: public GenericMapperNoExt<ContrTag, ContributionType, contributiontype_mapping>
{ public: using Base=GenericMapperNoExt<ContrTag,ContributionType,contributiontype_mapping>; using Base::str; };

#endif