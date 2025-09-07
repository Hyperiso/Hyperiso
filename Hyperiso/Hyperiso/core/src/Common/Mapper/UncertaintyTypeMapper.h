// #ifndef UNCERTAINTY_TYPE_MAPPER_H
// #define UNCERTAINTY_TYPE_MAPPER_H

// #include "GeneralEnum.h"
// #include "EnumMapper.h"

// class UncertaintyTypeMapper : public EnumMapperBase<UncertaintyType, UncertaintyTypeMapper> {
// public:
//     static const std::map<UncertaintyType, std::string>& mapping() {
//         static const std::map<UncertaintyType, std::string> m = {
//             {UncertaintyType::STAT, "Statistical"},
//             {UncertaintyType::SYST, "Systematics"},
//             {UncertaintyType::COMBINED, "Combined"}
//         };
//         return m;
//     }

//     static const std::map<std::string, UncertaintyType>& inverse_mapping() {
//         static const std::map<std::string, UncertaintyType> inv = invert_map(mapping());
//         return inv;
//     }

//     static const std::map<UncertaintyType, DataType>& data_type_mapping() {
//         static const std::map<UncertaintyType, DataType> m = {
//             {UncertaintyType::STAT, DataType::STD_STAT},
//             {UncertaintyType::SYST, DataType::STD_SYST},
//             {UncertaintyType::COMBINED, DataType::STD_COMBINED}
//         };
//         return m;
//     }

//     static DataType d_type(UncertaintyType u_type) {
//         return data_type_mapping().at(u_type);
//     }
// };

// #endif

#pragma once
#include "uncertaintytype_ids.hpp"
