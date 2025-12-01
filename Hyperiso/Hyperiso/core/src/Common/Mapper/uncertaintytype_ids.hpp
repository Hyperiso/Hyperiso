

#pragma once
#include "generic_mapper.hpp"
#include "Map.h"
#include "GeneralEnum.h"

struct UncTag {};
using UncertaintyTypeId = IdOf<UncTag>;

inline const std::map<UncertaintyType, std::string>& uncertaintytype_mapping() {
    static const std::map<UncertaintyType, std::string> m = {
        {UncertaintyType::STAT,     "STATISTICAL"},
        {UncertaintyType::SYST,     "SYSTEMATICS"},
        {UncertaintyType::COMBINED, "COMBINED"},
    };
    return m;
}

class UncertaintyTypeMapper
: public GenericMapperNoExt<UncTag, UncertaintyType, uncertaintytype_mapping>
{
public:
    using Base = GenericMapperNoExt<UncTag, UncertaintyType, uncertaintytype_mapping>;
    using Base::str;

    static const std::map<UncertaintyType, DataType>& data_type_mapping(){
        static const std::map<UncertaintyType, DataType> m = {
            {UncertaintyType::STAT,     DataType::STD_STAT},
            {UncertaintyType::SYST,     DataType::STD_SYST},
            {UncertaintyType::COMBINED, DataType::STD_COMBINED},
        };
        return m;
    }
    static DataType d_type(UncertaintyType t){ return data_type_mapping().at(t); }
};
