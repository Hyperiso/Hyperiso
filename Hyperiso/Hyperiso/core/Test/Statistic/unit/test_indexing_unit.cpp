#include "Include.h"
#include "Indexing.h"
#include "ObservableValue.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

static bool approx(double a, double b, double eps = 1e-12) {
    return std::fabs(a - b) <= eps;
}

int main() {
    std::cout << "== Indexing UNIT ==\n";

    {
        std::vector<int> ids{3, 1, 2};
        std::vector<double> vals{30.0, 10.0, 20.0};

        auto m = zip(ids, vals);
        assert(m.size() == 3);
        assert(approx(m.at(1), 10.0));
        assert(approx(m.at(2), 20.0));
        assert(approx(m.at(3), 30.0));

        auto u = unzip(m);
        assert((u.ids == std::vector<int>{1, 2, 3}));
        assert(u.vals.size() == 3);
        assert(approx(u.vals[0], 10.0));
        assert(approx(u.vals[1], 20.0));
        assert(approx(u.vals[2], 30.0));
    }

    {
        bool threw = false;
        try {
            (void)zip(std::vector<int>{1, 2}, std::vector<double>{1.0});
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    {
        std::vector<int> ids{10, 20};
        std::vector<std::vector<double>> vals{{1.0, 2.0}, {3.0, 4.0}};

        auto m = zip(ids, vals);
        assert(m.size() == 2);
        assert(approx(m.at(10).at(10), 1.0));
        assert(approx(m.at(10).at(20), 2.0));
        assert(approx(m.at(20).at(10), 3.0));
        assert(approx(m.at(20).at(20), 4.0));

        auto u = unzip(m);
        assert((u.ids == std::vector<int>{10, 20}));
        assert(u.vals.size() == 2);
        assert(approx(u.vals[0][0], 1.0));
        assert(approx(u.vals[0][1], 2.0));
        assert(approx(u.vals[1][0], 3.0));
        assert(approx(u.vals[1][1], 4.0));
    }

    {
        bool threw_empty = false;
        try {
            (void)zip(std::vector<int>{}, std::vector<std::vector<double>>{});
        } catch (const std::invalid_argument&) {
            threw_empty = true;
        }
        assert(threw_empty);

        bool threw_shape = false;
        try {
            (void)zip(std::vector<int>{1, 2}, std::vector<std::vector<double>>{{1.0, 2.0}});
        } catch (const std::invalid_argument&) {
            threw_shape = true;
        }
        assert(threw_shape);
    }

    {
        RealMatrix M(2, 2);
        M.at(0, 0) = 5.0;
        M.at(0, 1) = 6.0;
        M.at(1, 0) = 7.0;
        M.at(1, 1) = 8.0;

        auto m = zip(std::vector<int>{4, 9}, M);
        assert(approx(m.at(4).at(4), 5.0));
        assert(approx(m.at(4).at(9), 6.0));
        assert(approx(m.at(9).at(4), 7.0));
        assert(approx(m.at(9).at(9), 8.0));

        bool threw = false;
        try {
            (void)zip(std::vector<int>{1, 2, 3}, M);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        assert(threw);
    }

    {
        const ObservableId oid = ObservableMapper::to_id(Observables::BR_BS_MUMU);

        ObservableValue v0(oid, 1.25);
        ObservableValue v1(oid, 2.50);
        v1.bin = std::pair<double, double>{1.0, 6.0};

        std::map<ObservableId, std::vector<ObservableValue>> indexed;
        indexed.emplace(oid, std::vector<ObservableValue>{v0, v1});

        auto flat = flatten(indexed);
        assert(flat.ids.size() == 2);
        assert(flat.vals.size() == 2);
        assert(flat.ids[0].s == oid);
        assert(flat.ids[1].s == oid);
        assert(approx(flat.ids[0].p.first, 0.0));
        assert(approx(flat.ids[0].p.second, 0.0));
        assert(approx(flat.ids[1].p.first, 1.0));
        assert(approx(flat.ids[1].p.second, 6.0));
        assert(approx(flat.vals[0], 1.25));
        assert(approx(flat.vals[1], 2.50));
    }

    std::cout << "\nIndexing UNIT tests passed!\n";
    return 0;
}
