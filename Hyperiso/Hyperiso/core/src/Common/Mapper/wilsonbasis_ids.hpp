#pragma once
#include "generic_mapper.hpp"
#include "Map.h"

struct WilsonBasisTag {};
using WilsonBasisId = IdOf<WilsonBasisTag>;

class WilsonBasisMapper
: public GenericMapperNoExt<WilsonBasisTag, WilsonBasis, wilsonbasis_mapping>
{ public: using Base=GenericMapperNoExt<WilsonBasisTag,WilsonBasis,wilsonbasis_mapping>; using Base::str; };
