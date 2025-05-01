#ifndef QCD_ORDER_MAPPER_H
#define QCD_ORDER_MAPPER_H

#include "GeneralEnum.h"
#include "EnumMapper.h"

class OrderMapper : public EnumMapperBase<QCDOrder, OrderMapper> {
public:
    static const std::map<QCDOrder, std::string>& mapping() {
        static const std::map<QCDOrder, std::string> m = {
            {QCDOrder::NONE, "None"},
            {QCDOrder::LO, "LO"},
            {QCDOrder::NLO, "NLO"},
            {QCDOrder::NNLO, "NNLO"}
        };
        return m;
    }

    static const std::map<std::string, QCDOrder>& inverse_mapping() {
        static const std::map<std::string, QCDOrder> inv = invert_map(mapping());
        return inv;
    }
};

#endif