#pragma once
#include "generic_mapper.hpp"
#include "Map.h"

struct QCDOrderTag {};
using QCDOrderId = IdOf<QCDOrderTag>;

class OrderMapper
: public GenericMapperNoExt<QCDOrderTag, QCDOrder, order_mapping>
{
public:
    using Base = GenericMapperNoExt<QCDOrderTag, QCDOrder, order_mapping>;

    using Base::str;
    using Base::id_of;
    using Base::to_id;
    using Base::enum_of;
    using Base::list_all;
    using Base::init_builtins;

    static QCDOrder enum_elt(std::string_view s) {
        const auto key = normalize_key(s);
        for (const auto& [e, name] : order_mapping()) {
            if (normalize_key(name) == key) return e;
        }
        throw std::out_of_range("Unknown QCD order: " + std::string(s));
    }
};