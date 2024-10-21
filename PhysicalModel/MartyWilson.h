#pragma once
#include <complex>
#include "Wilsonv2.h"
#include "Interpolator.h"
#include "DataFrame.h"
#include "CSVReader.h"
#include "config.hpp"
#include <iostream>

class MartyWilson : public WilsonCoefficient {
public:
    MartyWilson(double Q_match, const std::string& coeff_name, const std::string& csv_path)
        : WilsonCoefficient(Q_match) {
        this->set_name(coeff_name);
        df = csv_reader.read_csv(csv_path);
        df.setIndex(df.getColumn<double>("Q_match").to_string_vec());
    }

    MartyWilson(double Q_match, const std::string& coeff_name)
        : WilsonCoefficient(Q_match) {
        std::string csv_path = project_root.data()+std::string("") + "DataBase/MartyWilson/SM_wilson.csv";
        this->set_name(coeff_name);
        df = csv_reader.read_csv(csv_path);
        df.setIndex(df.getColumn<double>("Q_match").to_string_vec());
    }

    std::complex<double> LO_calculation() override {
        auto closest_indices = find_closest_Q_matches(get_Q_match());

        double Q1 = df.iat<double>(closest_indices.first, "Q_match");
        double Q2 = df.iat<double>(closest_indices.second, "Q_match");

        std::complex<double> C1 = {
            df.iat<double>(closest_indices.first, this->get_name() + "_real"),
            df.iat<double>(closest_indices.first, this->get_name() + "_img")
        };

        std::complex<double> C2 = {
            df.iat<double>(closest_indices.second, this->get_name() + "_real"),
            df.iat<double>(closest_indices.second, this->get_name() + "_img")
        };

        return Interpolator::linearInterpolation(Q1, Q2, get_Q_match(), C1, C2);
    }

    std::complex<double> NLO_calculation() override {}
    std::complex<double> NNLO_calculation() override {}

private:
    CSVReader csv_reader;
    DataFrame df;

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

class BCoefficientGroupMarty : public BCoefficientGroup {

    BCoefficientGroupMarty(double Q_match) {
        this->insert(std::make_pair("C1", std::make_shared<MartyWilson>(Q_match, "C1"))); this->insert(std::make_pair("C2", std::make_shared<MartyWilson>(Q_match, "C2"))); this->insert(std::make_pair("C3", std::make_shared<MartyWilson>(Q_match, "C3")));
        this->insert(std::make_pair("C4", std::make_shared<MartyWilson>(Q_match, "C4")));  this->insert(std::make_pair("C5", std::make_shared<MartyWilson>(Q_match, "C5"))); this->insert(std::make_pair("C6", std::make_shared<MartyWilson>(Q_match, "C6"))); 
        this->insert(std::make_pair("C7", std::make_shared<MartyWilson>(Q_match, "C7")));  this->insert(std::make_pair("C8", std::make_shared<MartyWilson>(Q_match, "C8")));  this->insert(std::make_pair("C9", std::make_shared<MartyWilson>(Q_match, "C7"))); 
        this->insert(std::make_pair("C10", std::make_shared<MartyWilson>(Q_match, "C10")));
    }
};