#ifndef PARAMETER_TYPE_IDS_H
#define PARAMETER_TYPE_IDS_H

#include "generic_mapper.h"
#include "Map.h"

/**
 * @file ParameterTypeMapper.h
 * @brief Mapper for ParameterType <-> string identifiers.
 */
struct ParamTypeTag {};
using ParameterTypeId = IdOf<ParamTypeTag>;

class ParameterTypeMapper
: public GenericMapperNoExt<ParamTypeTag, ParameterType, parametertype_mapping>
{ public: using Base=GenericMapperNoExt<ParamTypeTag,ParameterType,parametertype_mapping>; using Base::str; };

#endif