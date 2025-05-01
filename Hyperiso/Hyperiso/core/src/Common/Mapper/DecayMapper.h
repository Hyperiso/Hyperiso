#ifndef DECAY_MAPPER_H
#define DECAY_MAPPER_H


#include "GeneralEnum.h"
#include "EnumMapper.h"
#include "Logger.h"

class DecayMapper : public EnumMapperBase<Decays, DecayMapper> {
public:
    static const std::map<Decays, std::string>& mapping();
    static const std::map<std::string, Decays>& inverse_mapping();

    static std::vector<Observables> get_observables(Decays decay);

    static Decays get_decay(Observables obs);
};

#endif