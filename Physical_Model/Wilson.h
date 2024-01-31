#if !defined(HYPERISO_WILSON_H)
#define HYPERISO_WILSON_H

#include <complex>
#include <vector>
#include <string>
#include <memory>
#include "./Math/Math.h"
#include "QCDParameters.h"
#include "Core/Parameters.h"
#include "Wilson_parameters.h"
#define MAX_ORDER 2
#define NW 18

typedef std::complex<double> complex_t; 
typedef std::vector<std::vector<complex_t>> WilsonSet;

enum class WilsonCoefficient {
    C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, CQ1, CQ2, CP7, CP8, CP9, CP10, CPQ1, CPQ2
};

class WilsonManager{
private:
    static WilsonManager* instance;
    QCDParameters run;
    WilsonSet C_match;  // C1...10, CQ1, CQ2, C'7...10, CQ'1, CQ'2 for each order
    WilsonSet C;
    const double Q_matching;
    double Q;
    std::unique_ptr<InitializationStrategy> strategy;

    inline explicit WilsonManager(double mu_match, InitializationStrategy* strategy): Q_matching(mu_match) {
        WilsonInitializer wi{mu_match, strategy};
        wi.init(C_match);
    }

public:
    WilsonManager(WilsonManager&) = delete;
    void operator=(const WilsonManager&) = delete;
    static WilsonManager* GetInstance(double mu_match=0, InitializationStrategy* strategy = nullptr) {
        if (!WilsonManager::instance) {
            if (mu_match == 0.0) {
                // Log an error
            } else {
                WilsonManager::instance = new WilsonManager(mu_match, strategy);
            }
        }
        return WilsonManager::instance;
    }

    complex_t get(WilsonCoefficient wc, int order) const;    // Returns C_id at a given order
    void setScale(double mu) {
        strategy->set_base1(C, C_match, Q, Q_matching, run);  // Appel à set_base1
    }               // Computes the C's at scale mu using RGEs 
};  

class WilsonInitializer {
private:
    Parameters* sm = Parameters::GetInstance();
    std::unique_ptr<InitializationStrategy> strategy;
    double scale;
    QCDParameters run;
    void scanLHA(std::string_view file); 
    // void initSM_LO();
    // void initSM_NLO();
    // void initSM_NNLO();
    // void initTHDM();
    // void initSUSY();

public:
    WilsonInitializer(double mu_match, InitializationStrategy* strat) :
    scale(mu_match), strategy(strat) {}

    void init(WilsonSet& C_match) {
        strategy->init(sm, scale, C_match, run);
    }
    
};

class InitializationStrategy {

    

public:
    virtual void init(Parameters* sm, double scale, WilsonSet& C_match, QCDParameters& run) = 0;
    virtual void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match, QCDParameters& run) = 0;
    virtual void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match, QCDParameters& run) =0;
    virtual ~InitializationStrategy() {}
};

class SM_LO_Strategy : public InitializationStrategy {
public:
    void init(Parameters* sm, double scale, WilsonSet& C_match, QCDParameters& run) override;
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match, QCDParameters& run) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match, QCDParameters& run) override;
};


class SM_NLO_Strategy : public InitializationStrategy {
public:
    void init(Parameters* sm, double scale, WilsonSet& C_match, QCDParameters& run) override;
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match, QCDParameters& run) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match, QCDParameters& run) override;
};

class SM_NNLO_Strategy : public InitializationStrategy {
public:
    void init(Parameters* sm, double scale, WilsonSet& C_match, QCDParameters& run) override;
    void set_base1(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match, QCDParameters& run) override;
    void set_base2(WilsonSet& C, WilsonSet& C_match, double Q, const double Q_match, QCDParameters& run) override;
};
#endif // HYPERISO_WILSON_H
