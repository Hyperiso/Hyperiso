#pragma once
#include "generic_mapper.hpp"
#include "Map.h"

struct ModelTag {};
using ModelId = IdOf<ModelTag>;

class ModelMapper
: public GenericMapperNoExt<ModelTag, Model, model_mapping>
{ public: using Base=GenericMapperNoExt<ModelTag,Model,model_mapping>; using Base::str; };
