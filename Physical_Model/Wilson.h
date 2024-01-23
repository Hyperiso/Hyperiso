#if !defined(HYPERISO_WILSON_H)
#define HYPERISO_WILSON_H

#include <complex>
#include <vector>

#define MAX_ORDER 2
#define NW 18

typedef std::complex<double> complex_t; 
typedef std::vector<std::vector<complex_t>> WilsonSet;

class WilsonManager{
private:
    WilsonSet C_match;  // C1...10, CQ1, CQ2, C'7...10, CQ'1, CQ'2 for each order
    WilsonSet C;
    const double matching_scale;
    double scale;

public:
    inline explicit WilsonManager(double mu_match): matching_scale(mu_match) {
        WilsonInitializer wi{C_match, mu_match};
        wi.init();
    }
    complex_t get(size_t id, int order);    // Returns C_id at a given order
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
