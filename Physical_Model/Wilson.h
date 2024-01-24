#if !defined(HYPERISO_WILSON_H)
#define HYPERISO_WILSON_H

#include <complex>
#include <vector>

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

    WilsonSet C_match;  // C1...10, CQ1, CQ2, C'7...10, CQ'1, CQ'2 for each order
    WilsonSet C;
    const double matching_scale;
    double scale;

    inline explicit WilsonManager(double mu_match): matching_scale(mu_match) {
        WilsonInitializer wi{C_match, mu_match};
        wi.init();
    }

public:
    WilsonManager(WilsonManager&) = delete;
    void operator=(const WilsonManager&) = delete;
    static WilsonManager* GetInstance(double mu_match=0);

    complex_t get(WilsonCoefficient wc, int order) const;    // Returns C_id at a given order
    void setScale(double mu);               // Computes the C's at scale mu using RGEs  
};  

class WilsonInitializer {
private:
    double scale;
    WilsonSet C;

    void scanLHA(std::string_view file); 
    void initSM_LO();
    void initSM_NLO();
    void initSM_NNLO();
    void initTHDM();
    void initSUSY();

public:
    WilsonInitializer(WilsonSet& C_empty, double mu_match) : C(C_empty), scale(mu_match) {}
    void init();
};

#endif // HYPERISO_WILSON_H
