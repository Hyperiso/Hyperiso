#pragma once
#include "generic_mapper.hpp"
#include "Map.h"

struct ParamTypeTag {};
using ParameterTypeId = IdOf<ParamTypeTag>;

class ParameterTypeMapper
: public GenericMapperNoExt<ParamTypeTag, ParameterType, parametertype_mapping>
{ public: using Base=GenericMapperNoExt<ParamTypeTag,ParameterType,parametertype_mapping>; using Base::str; };
