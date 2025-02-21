#pragma once
#include <complex>
#include "Wilson.h"
#include "Interpolator.h"
#include "DataFrame.h"
#include "CSVReader.h"
#include "MartyInterface.h"
#include "config.hpp"
#include <iostream>
#include <math.h>

class MartyWilson : public WilsonCoefficient {
public:
    MartyWilson(double Q_match, const std::string& coeff_name, const std::string& csv_path)
        : WilsonCoefficient(Q_match) {
        this->csv_path = csv_path;
        this->set_name(coeff_name);
        df = csv_reader.read_csv(csv_path);
        df.setIndex(df.getColumn<double>("Q_match").to_string_vec());
    }

    MartyWilson(double Q_match, const std::string& coeff_name)
        : WilsonCoefficient(Q_match) {
        this->set_name(coeff_name);
        df = csv_reader.read_csv(this->csv_path);
        df.setIndex(df.getColumn<double>("Q_match").to_string_vec());
    }

    MartyWilson(const std::string& coeff_name)
        : WilsonCoefficient() {
        this->set_name(coeff_name);
        df = csv_reader.read_csv(this->csv_path);
        df.setIndex(df.getColumn<double>("Q_match").to_string_vec());
    }

    std::string get_model() {
        return this->model;
    }
    void set_model(std::string model) {
        this->model = model;
    }

    std::complex<double> LO_calculation() override {
        double epsi = 1e-4;
        for (size_t i = 0; i < df.getRowCount(); ++i) {
            double Q_match = df.iat<double>(i, "Q_match");
            if (fabs(Q_match-this->get_Q_match()) < epsi) {
                std::cout << this->get_name() << " waw" << std::endl;
                for (auto& _ : this->df.getColumnNames()) {
                    if (this->get_name()+"_real" == _) {
                        if (isnan(df.iat<double>(i, this->get_name()+"_real")) && isnan(df.iat<double>(i, this->get_name()+"_img"))) {
                            break;
                        }
                        std::cout << df.iat<double>(i, this->get_name()+"_real") << " BUTE" << std::endl;
                        this->set_CoefficientMatchingValue("LO", {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")});
                        return {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")};
                    }
                } 
            }
        }
        MartyInterface MartyInterface;
        MartyInterface.calculate(this->get_name(), this->get_model(), this->get_Q_match());
        df = csv_reader.read_csv(this->csv_path);
        df.setIndex(df.getColumn<double>("Q_match").to_string_vec());


        for (size_t i = 0; i < df.getRowCount(); ++i) {
            double Q_match = df.iat<double>(i, "Q_match");
            if (fabs(Q_match-this->get_Q_match()) < epsi) {
                this->set_CoefficientMatchingValue("LO", {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")});
                return {df.iat<double>(i, this->get_name()+"_real"), df.iat<double>(i, this->get_name()+"_img")};
            }
        }

        return {0., 0.};
    }

    std::complex<double> NLO_calculation() override {return {0., 0.};}
    std::complex<double> NNLO_calculation() override {return {0., 0.};}

private:
    CSVReader csv_reader;
    DataFrame df;
    std::string model{"SM"};
    std::string csv_path{project_root.data() + std::string("") + "../MartyWilson/SM_wilson.csv"};
    std::pair<size_t, size_t> find_closest_Q_matches(double target_Q_match) {
        
        size_t closest_below = 0, closest_above = 0;
        double min_diff_below = std::numeric_limits<double>::max();
        double min_diff_above = std::numeric_limits<double>::max();

        for (size_t i = 0; i < df.getRowCount(); ++i) {
            double Q_match = df.iat<double>(i, "Q_match");
            double diff = Q_match - target_Q_match;
            std::cout << "Q_match : " << Q_match << std::endl;
            if (diff <= 0 && std::abs(diff) < min_diff_below) {
                closest_below = i;
                min_diff_below = std::abs(diff);
            } else if (diff > 0 && diff < min_diff_above) {
                closest_above = i;
                min_diff_above = diff;
            }
        }
        std::cout << closest_above << " " << closest_below << std::endl;
        return {closest_below, closest_above};
    }
};

// class BCoefficientGroupMarty : public BCoefficientGroup {
// public:
//     BCoefficientGroupMarty(double Q_match) { this->clear();
//         this->insert(std::make_pair("C1", std::make_shared<MartyWilson>(Q_match, "C1"))); this->insert(std::make_pair("C2", std::make_shared<MartyWilson>(Q_match, "C2"))); this->insert(std::make_pair("C3", std::make_shared<MartyWilson>(Q_match, "C3")));
//         this->insert(std::make_pair("C4", std::make_shared<MartyWilson>(Q_match, "C4")));  this->insert(std::make_pair("C5", std::make_shared<MartyWilson>(Q_match, "C5"))); this->insert(std::make_pair("C6", std::make_shared<MartyWilson>(Q_match, "C6"))); 
//         this->insert(std::make_pair("C7", std::make_shared<MartyWilson>(Q_match, "C7")));  this->insert(std::make_pair("C8", std::make_shared<MartyWilson>(Q_match, "C8")));  this->insert(std::make_pair("C9", std::make_shared<MartyWilson>(Q_match, "C7"))); 
//         this->insert(std::make_pair("C10", std::make_shared<MartyWilson>(Q_match, "C10")));
//     }

//     BCoefficientGroupMarty() { this->clear();
//         this->insert(std::make_pair("C1", std::make_shared<MartyWilson>(81, "C1"))); this->insert(std::make_pair("C2", std::make_shared<MartyWilson>(81, "C2"))); this->insert(std::make_pair("C3", std::make_shared<MartyWilson>(81, "C3")));
//         this->insert(std::make_pair("C4", std::make_shared<MartyWilson>(81, "C4")));  this->insert(std::make_pair("C5", std::make_shared<MartyWilson>(81, "C5"))); this->insert(std::make_pair("C6", std::make_shared<MartyWilson>(81, "C6"))); 
//         this->insert(std::make_pair("C7", std::make_shared<MartyWilson>(81, "C7")));  this->insert(std::make_pair("C8", std::make_shared<MartyWilson>(81, "C8")));  this->insert(std::make_pair("C9", std::make_shared<MartyWilson>(81, "C7"))); 
//         this->insert(std::make_pair("C10", std::make_shared<MartyWilson>(81, "C10")));
//     }
// };

// class BPrimeCoefficientGroupMarty : public BPrimeCoefficientGroup {
// public:
//     BPrimeCoefficientGroupMarty(double Q_match) { this->clear();
//         this->insert(std::make_pair("C1", std::make_shared<MartyWilson>(Q_match, "C1"))); this->insert(std::make_pair("C2", std::make_shared<MartyWilson>(Q_match, "C2"))); this->insert(std::make_pair("C3", std::make_shared<MartyWilson>(Q_match, "C3")));
//         this->insert(std::make_pair("C4", std::make_shared<MartyWilson>(Q_match, "C4")));  this->insert(std::make_pair("C5", std::make_shared<MartyWilson>(Q_match, "C5"))); this->insert(std::make_pair("C6", std::make_shared<MartyWilson>(Q_match, "C6"))); 
//         this->insert(std::make_pair("C7", std::make_shared<MartyWilson>(Q_match, "C7")));  this->insert(std::make_pair("C8", std::make_shared<MartyWilson>(Q_match, "C8")));  this->insert(std::make_pair("C9", std::make_shared<MartyWilson>(Q_match, "C7"))); 
//         this->insert(std::make_pair("C10", std::make_shared<MartyWilson>(Q_match, "C10")));
//     }

//     BPrimeCoefficientGroupMarty() { this->clear();
//         this->insert(std::make_pair("C1", std::make_shared<MartyWilson>(81, "C1"))); this->insert(std::make_pair("C2", std::make_shared<MartyWilson>(81, "C2"))); this->insert(std::make_pair("C3", std::make_shared<MartyWilson>(81, "C3")));
//         this->insert(std::make_pair("C4", std::make_shared<MartyWilson>(81, "C4")));  this->insert(std::make_pair("C5", std::make_shared<MartyWilson>(81, "C5"))); this->insert(std::make_pair("C6", std::make_shared<MartyWilson>(81, "C6"))); 
//         this->insert(std::make_pair("C7", std::make_shared<MartyWilson>(81, "C7")));  this->insert(std::make_pair("C8", std::make_shared<MartyWilson>(81, "C8")));  this->insert(std::make_pair("C9", std::make_shared<MartyWilson>(81, "C7"))); 
//         this->insert(std::make_pair("C10", std::make_shared<MartyWilson>(81, "C10")));
//     }
// };

// class BScalarCoefficientGroupMarty : public BScalarCoefficientGroup {
// public:
//     BScalarCoefficientGroupMarty(double Q_match) { this->clear();
//         this->insert(std::make_pair("C1", std::make_shared<MartyWilson>(Q_match, "C1"))); this->insert(std::make_pair("C2", std::make_shared<MartyWilson>(Q_match, "C2"))); this->insert(std::make_pair("C3", std::make_shared<MartyWilson>(Q_match, "C3")));
//         this->insert(std::make_pair("C4", std::make_shared<MartyWilson>(Q_match, "C4")));  this->insert(std::make_pair("C5", std::make_shared<MartyWilson>(Q_match, "C5"))); this->insert(std::make_pair("C6", std::make_shared<MartyWilson>(Q_match, "C6"))); 
//         this->insert(std::make_pair("C7", std::make_shared<MartyWilson>(Q_match, "C7")));  this->insert(std::make_pair("C8", std::make_shared<MartyWilson>(Q_match, "C8")));  this->insert(std::make_pair("C9", std::make_shared<MartyWilson>(Q_match, "C7"))); 
//         this->insert(std::make_pair("C10", std::make_shared<MartyWilson>(Q_match, "C10")));
//     }

//     BScalarCoefficientGroupMarty() { this->clear();
//         this->insert(std::make_pair("C1", std::make_shared<MartyWilson>(81, "C1"))); this->insert(std::make_pair("C2", std::make_shared<MartyWilson>(81, "C2"))); this->insert(std::make_pair("C3", std::make_shared<MartyWilson>(81, "C3")));
//         this->insert(std::make_pair("C4", std::make_shared<MartyWilson>(81, "C4")));  this->insert(std::make_pair("C5", std::make_shared<MartyWilson>(81, "C5"))); this->insert(std::make_pair("C6", std::make_shared<MartyWilson>(81, "C6"))); 
//         this->insert(std::make_pair("C7", std::make_shared<MartyWilson>(81, "C7")));  this->insert(std::make_pair("C8", std::make_shared<MartyWilson>(81, "C8")));  this->insert(std::make_pair("C9", std::make_shared<MartyWilson>(81, "C7"))); 
//         this->insert(std::make_pair("C10", std::make_shared<MartyWilson>(81, "C10")));
//     }
// };