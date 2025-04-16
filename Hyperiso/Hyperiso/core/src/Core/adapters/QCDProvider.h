#ifndef __QCDADAPTER_H__
#define __QCDADAPTER_H__

#include "IDataProvider.h"
#include "IQCDProvider.h"
#include "QCDHelper.h"

struct AlphasConfig : public AbstractConfig {
    double scale;
    MassType m_b_type;
    MassType m_t_type;

    AlphasConfig(double scale, MassType m_b_type, MassType m_t_type) : scale(scale), m_b_type(m_b_type), m_t_type(m_t_type) {}
};

struct MassConfig : public AlphasConfig {
    int pdg_id;

    MassConfig(int pdg_id, double scale, MassType m_b_type, MassType m_t_type) : AlphasConfig(scale, m_b_type, m_t_type), pdg_id(pdg_id) {}
};

class QCDProvider : public IDataProvider<QCDProvider>, public IQCDProvider {
public:
    double operator()(AlphasConfig);
    double operator()(MassConfig);
    QCDConstants* get_constants() override;
};


#endif // __QCDADAPTER_H__
