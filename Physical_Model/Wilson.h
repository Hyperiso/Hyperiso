#if !defined(HYPERISO_WILSON_H)
#define HYPERISO_WILSON_H

#include <complex>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include "../Math/Math.h"
#include "../Core/Parameters.h"
#include "../Core/Logger.h"
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
    virtual void init(Parameters* sm, double scale, WilsonSet& C_match) = 0;
    virtual void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) = 0;
    virtual void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) =0;
    virtual ~InitializationStrategy() {}
};

class SM_LO_Strategy : public InitializationStrategy {
public:
    void init(Parameters* sm, double scale, WilsonSet& C_match) override;
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
};


class SM_NLO_Strategy : public SM_LO_Strategy {
public:
    void init(Parameters* sm, double scale, WilsonSet& C_match) override;
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match) override;
};

class SM_NNLO_Strategy : public SM_NLO_Strategy {
public:
    void init(Parameters* sm, double scale, WilsonSet& C_match) override;
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
        strategy->init(sm, scale, C_match);
    }
    
};




class WilsonManager{
private:
    static std::map<std::string, WilsonManager*> instances;
    WilsonSet C_match;  // C1...10, CQ1, CQ2, C'7...10, CQ'1, CQ'2 for each order
    WilsonSet C;
    const double Q_matching;
    double Q;
    std::shared_ptr<InitializationStrategy> strategy;

    inline explicit WilsonManager(double mu_match, std::shared_ptr<InitializationStrategy> strat): Q_matching(mu_match), strategy(strat) {
        
        std::cout <<"Creation of WilsonManager" << std::endl;
        WilsonInitializer wi{mu_match, strategy};
        std::cout <<"Creation of WilsonInitializer Done" << std::endl;
        wi.init(C_match);
        // std::cout << "First wilson Coefficient " << C_match[0][static_cast<size_t>(WilsonCoefficient::C7)] << std::endl;
    }

public:
    WilsonManager(WilsonManager&) = delete;
    void operator=(const WilsonManager&) = delete;
    static WilsonManager* GetInstance(const std::string& strategyName, double mu_match, std::shared_ptr<InitializationStrategy> strategy) {
        // Vérifier si l'instance existe déjà pour cette stratégie
        auto it = instances.find(strategyName);
        if (it == instances.end()) {
            // Vérification de mu_match
            if (mu_match == 0.0) {
                std::cerr << "Erreur: mu_match ne doit pas être 0 lors de la première création de WilsonManager." << std::endl;
                std::exit(EXIT_FAILURE);
            }
            // Créer une nouvelle instance
            WilsonManager* newInstance = new WilsonManager(mu_match, strategy);
            instances[strategyName] = newInstance;
            return newInstance;
        }
        // Retourner l'instance existante
        return it->second;
    }

    static void Cleanup() {
        for (auto& pair : instances) {
            delete pair.second;
        }
        instances.clear();
    }

    complex_t get(WilsonCoefficient wc, int order) const;    // Returns C_id at a given order
    complex_t get_matchs(WilsonCoefficient wc_match, int order) const;
    void setScale(double mu) {
        Q= mu;
        strategy->set_base1(C, C_match, Q, Q_matching);  // Appel à set_base1
    }               // Computes the C's at scale mu using RGEs 
};  



#endif // HYPERISO_WILSON_H
