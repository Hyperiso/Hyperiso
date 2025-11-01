#include "WilsonGroup.h"

class KCoefficientGroup : public CoefficientGroup {
public:
    KCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm=false);

    void set_gen(int new_gen) {}
    std::shared_ptr<CoefficientGroup> clone() const override;
    
    std::shared_ptr<CoefficientGroup> get_sm_group() override { return std::make_shared<KCoefficientGroup>(adapters, true); }
    static std::unordered_map<WCoef, scalar_t> base_1_LO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    // static std::unordered_map<WCoef, scalar_t> base_2_LO_calculation   (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    // static std::unordered_map<WCoef, scalar_t> base_1_NLO_calculation  (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    // static std::unordered_map<WCoef, scalar_t> base_2_NLO_calculation  (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    // static std::unordered_map<WCoef, scalar_t> base_1_NNLO_calculation (const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching, const std::unordered_map<std::string, std::shared_ptr<Block>>& src);
    
    private:

    // void init_running_parameter_blocks();

    void init_sources() {}
    void add_wilson_coefficients(bool force_sm=false) {}
};