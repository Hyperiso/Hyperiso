#include "Logger.h"
#include <iostream>
#include "config.hpp"
// #include "hpl.h"
#include "Math.h"
#include "BXsllDecay.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    // auto wb = std::shared_ptr<ObsWilsonBuilder>(nullptr);
    // BXsllDecay d {QCDOrder::LO, 0., 0., wb};

    // constexpr size_t N = 10000;
    // std::ofstream of;
    // of.open("Bsll_func.csv");

    // auto write_line = [&of, &d] (scalar_t s) {
    //     double sd = real(s);
    //     of << sd << ","
    //        << real(log(s)) << ","
    //        << real(log(1.-s)) << ","
    //        << real(sLi2(s)) << ","
    //        << real(sLi3(s)) << ","
    //        << real(sLi3(1.-s)) << ","
    //        << real(sLi4(s)) << ","
    //        << real(sqrt(s)) << ","
    //        << real(sqrt(4.-s)) << ","
    //        << real(asin(sqrt(s)/2.)) << ","
    //        << real(sCl2(2.*asin(sqrt(sd)/2.))) << ","
    //        << real(sCl3(2.*asin(sqrt(sd)/2.))) << ","
    //        << real(d.f_17(s, 0, 0.0625)) << "," << imag(d.f_17(s, 0, 0.0625)) << ","
    //        << real(d.f_19(s, 0, 0.0625)) << "," << imag(d.f_19(s, 0, 0.0625)) << ","
    //        << real(d.f_27(s, 0, 0.0625)) << "," << imag(d.f_27(s, 0, 0.0625)) << ","
    //        << real(d.f_29(s, 0, 0.0625)) << "," << imag(d.f_29(s, 0, 0.0625)) 
    //        << "\n";
    // };

    // of << "s_hat,Ls,Lsb,Li2s,Li3s,Li3sb,Li4s,sqrts,sqrt4s,as,cl2,cl3,f17r,f17i,f19r,f19i,f27r,f27i,f29r,f29i\n";
    // // write_line(1e-6);
    // for (std::size_t i = 1; i < N - 1; ++i) {
    //     double s = static_cast<double>(i) / static_cast<double>(N - 1);
    //     write_line(s);
    // }
    // // write_line(1 - 1e-6);
    // of.close();

    return 0;
}