#ifndef __MARTYADAPTER_H__
#define __MARTYADAPTER_H__

#include "IMonitor.h"
#include "IPathProvider.h"

enum class MartyPath {
    MODEL_FILE,
    TEMPLATE_DIR,
    PARAM_MAPPING_DIR,
};

class MartyAdapter : public IMonitor<InternalFlag>, IPathProvider<MartyPath> {
private:
    /* data */
public:
    fs::path get_path(MartyPath path_name) override;
    bool check_flag(InternalFlag flag) override;
    std::string get_marty_model_name();
};


#endif // __MARTYADAPTER_H__
