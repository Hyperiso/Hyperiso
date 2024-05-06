#if !defined(HYPERISO_WILSON_H)
#define HYPERISO_WILSON_H

#include <complex>
#include <vector>
#include <string>
#include <memory>
#include "Math.h"
#include "Parameters.h"
#include "Logger.h"
#include "Wilson_parameters.h"
#define MAX_ORDER 2
#define NW 18

typedef std::complex<double> complex_t; 
typedef std::vector<std::vector<complex_t>> WilsonSet;

enum class WilsonCoefficient {
    C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, CQ1, CQ2, CP7, CP8, CP9, CP10, CPQ1, CPQ2
};



class InitializationStrategy {

    

public:

    /**
     * @brief Initialize the Wilson coefficients at scale scale (usually 81 GeV).
     * 
     * @param sm Pointer to the Parameters object.
     * @param scale Matching scale.
     * @param C_match Wilson coefficients set at matching scale.
     */
    virtual void init(double scale, WilsonSet& C_match) = 0;

    /**
     * @brief Initialize the prime Wilson coefficients at scale scale (usually 81 GeV).
     * 
     * @param scale_W Matching scale.
     * @param scale Scale of interest.
     * @param gen number of generation.
     * @param C_match Wilson coefficient set at matching scale.
     */
    virtual void init_prime(double scale_W,double scale,int gen, WilsonSet& C_match) = 0;

    /**
     * @brief Initialize the scalar Wilson coefficients at scale scale (usually 81 GeV).
     * 
     * @param scale_W Matching scale.
     * @param scale Scale of interest.
     * @param gen number of generation.
     * @param C_match Wilson coefficient set at matching scale.
     */
    virtual void init_scalar(double scale_W,double scale,int gen, WilsonSet& C_match) = 0;

    /**
     * @brief Run the Wilson coefficient from the scale Q_match to the scale Q, in the conventional basis.
     * 
     * @param C Wilson coefficients at scale Q.
     * @param C_match Wilson coefficients at scale Q_match.
     * @param Q Scale.
     * @param Q_match Matching scale.
     */
    virtual void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) = 0;

    /**
     * @brief Run the Wilson coefficient from the scale Q_match to the scale Q, in the traditional basis.
     * 
     * @param C Wilson coefficients at scale Q.
     * @param C_match Wilson coefficients at scale Q_match.
     * @param Q Scale.
     * @param Q_match Matching scale.
     */
    virtual void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) =0;
    virtual ~InitializationStrategy() {}
};

class SM_LO_Strategy : public InitializationStrategy {
public:
    void init(double scale, WilsonSet& C_match) override;
    void init_prime(double scale_W,double scale,int gen, WilsonSet& C_match);
    void init_scalar(double scale_W,double scale,int gen, WilsonSet& C_match);
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
};


class SM_NLO_Strategy : public SM_LO_Strategy {
public:
    void init(double scale, WilsonSet& C_match) override;
    void init_prime(double scale_W,double scale,int gen, WilsonSet& C_match) {}
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
};

class SM_NNLO_Strategy : public SM_NLO_Strategy {
public:
    void init(double scale, WilsonSet& C_match) override;
    void init_prime(double scale_W,double scale,int gen, WilsonSet& C_match) {}
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override {}
};


class WilsonInitializer {
private:
    Parameters* sm = Parameters::GetInstance();
    std::shared_ptr<InitializationStrategy> strategy = std::make_shared<SM_LO_Strategy>();;
    double scale;


public:
    WilsonInitializer(double mu_match, std::shared_ptr<InitializationStrategy> strat) :
    scale(mu_match), strategy(strat){
        std::cout <<"Creation of WilsonManager" << std::endl;
    }

    void init(WilsonSet& C_match) {
        strategy->init(scale, C_match);
    }
    
};




class WilsonManager{
private:
    static std::map<std::string, WilsonManager*> instances;

    /**
     * @brief Matching Wilson coefficients.
     * 
     * Stores the matching Wilson coefficients for each order at scale Q_matching.
     */
    WilsonSet C_match;  // C1...10, CQ1, CQ2, C'7...10, CQ'1, CQ'2 for each order

    /**
     * @brief Wilson coefficients.
     * 
     * Stores the Wilson coefficients for each order at scale Q.
     */
    WilsonSet C;

    /**
     * @brief Matching scale.
     * 
     * The scale at which Wilson coefficients are matched.
     */
    const double Q_matching;

    /**
     * @brief Current scale.
     * 
     * The scale at which Wilson coefficients are computed.
     */
    double Q;
    std::shared_ptr<InitializationStrategy> strategy;

    /**
     * @brief Constructor for WilsonManager.
     * 
     * @param mu_match Matching scale.
     * @param strat Initialization strategy.
     */
    inline explicit WilsonManager(double mu_match, std::shared_ptr<InitializationStrategy> strat): Q_matching(mu_match), strategy(strat) {
        
        std::cout <<"Creation of WilsonManager" << std::endl;
        WilsonInitializer wi{mu_match, strategy};
        std::cout <<"Creation of WilsonInitializer Done" << std::endl;
        wi.init(C_match);
        // std::cout << "First wilson Coefficient " << C_match[0][static_cast<size_t>(WilsonCoefficient::C7)] << std::endl;
    }

public:
    /**
     * @brief Deleted copy constructor.
     */
    WilsonManager(WilsonManager&) = delete;

    /**
     * @brief Deleted copy assignment operator.
     */
    void operator=(const WilsonManager&) = delete;

    /**
     * @brief Get the instance of WilsonManager.
     * 
     * Creates a new instance if it does not exist for the given strategy, otherwise returns the existing instance.
     * 
     * @param strategyName Name of the strategy.
     * @param mu_match Matching scale.
     * @param strategy Initialization strategy.
     * @return Pointer to the WilsonManager instance.
     */
    static WilsonManager* GetInstance(const std::string& strategyName, double mu_match, std::shared_ptr<InitializationStrategy> strategy) {

        auto it = instances.find(strategyName);
        if (it == instances.end()) {
            if (mu_match == 0.0) {
                std::cerr << "Error: mu_match should not be 0 during the first creation of WilsonManager." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            WilsonManager* newInstance = new WilsonManager(mu_match, strategy);
            instances[strategyName] = newInstance;
            return newInstance;
        }
        return it->second;
    }

    /**
     * @brief Clean up the instances of WilsonManager.
     * 
     * Deletes all instances of WilsonManager and clears the map.
     */
    static void Cleanup() {
        for (auto& pair : instances) {
            delete pair.second;
        }
        instances.clear();
    }

    /**
     * @brief Get the Wilson coefficient at a given order.
     * 
     * @param wc Wilson coefficient.
     * @param order Order of the coefficient.
     * @return Wilson coefficient value.
     */
    complex_t get(WilsonCoefficient wc, int order) const;    // Returns C_id at a given order

    /**
     * @brief Get the Wilson coefficient up to a given order.
     * 
     * @param wc Wilson coefficient.
     * @param order Maximal order of the coefficient.
     * @return Wilson coefficient value.
     */
    complex_t get_full(WilsonCoefficient wc, int order) const;    // Returns C_id at a given order

    /**
     * @brief Get the matched Wilson coefficient at a given order.
     * 
     * @param wc_match Matched Wilson coefficient.
     * @param order Order of the coefficient.
     * @return Matched Wilson coefficient value.
     */
    complex_t get_matchs(WilsonCoefficient wc_match, int order) const;

    /**
     * @brief Set the scale for computing Wilson coefficients.
     * 
     * Computes the Wilson coefficients at the specified scale using the Renormalization Group Equations (RGEs).
     * 
     * @param mu Scale value.
     */
    void setScale(double mu) {
        Q= mu;
        strategy->set_base1(C, C_match, Q, Q_matching);
    }

};  



#endif // HYPERISO_WILSON_H
