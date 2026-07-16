#ifndef MODEL_IDS_H
#define MODEL_IDS_H

#include "generic_mapper.h"
#include "Map.h"

/**
 * @file ModelMapper.h
 * @brief Mapper for Model <-> string identifiers.
 */
struct ModelTag {};
using ModelId = IdOf<ModelTag>;

class ModelMapper
: public GenericMapperNoExt<ModelTag, Model, model_mapping>
{ public: using Base=GenericMapperNoExt<ModelTag,Model,model_mapping>; using Base::str; };

#endif