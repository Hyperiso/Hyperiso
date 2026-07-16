#include "Include.h"
#include "Indexing.h"
#include "ObservableValue.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

static bool approx(double a, double b, double eps = 1e-12) {
    return std::fabs(a - b) <= eps;
}

int main() {
    std::cout << "== Indexing INTEGRATION ==\n";

    {
        std::vector<std::size_t> ids{0, 2, 5};
        std::vector<double> vals{1.0, 4.0, 25.0};

        auto indexed = zip(ids, vals);
        auto dense = unzip(indexed);

        assert(dense.ids.size() == indexed.size());
        assert(dense.vals.size() == indexed.size());
        assert(dense.ids[0] == 0);
        assert(dense.ids[1] == 2);
        assert(dense.ids[2] == 5);
        assert(approx(dense.vals[0], 1.0));
        assert(approx(dense.vals[1], 4.0));
        assert(approx(dense.vals[2], 25.0));
    }

    {
        std::vector<int> ids{1, 3, 4};
        RealMatrix M(3, 3);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                M.at(i, j) = 10.0 * static_cast<double>(i + 1) + static_cast<double>(j + 1);
            }
        }

        auto indexed = zip(ids, M);
        auto dense = unzip(indexed);

        assert((dense.ids == ids));
        assert(dense.vals.size() == 3);
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                assert(approx(dense.vals[i][j], M.at(i, j)));
            }
        }
    }

    {
        const ObservableId obs = ObservableMapper::to_id(Observables::BR_BS_MUMU);

        ObservableValue v0(obs, 3.0);
        ObservableValue v1(obs, 4.0);
        v1.bin = std::pair<double, double>{0.1, 1.0};
        ObservableValue v2(obs, 5.0);
        v2.bin = std::pair<double, double>{1.0, 2.0};

        std::map<ObservableId, std::vector<ObservableValue>> indexed;
        indexed[obs] = {v0, v1, v2};

        auto flat = flatten(indexed);
        assert(flat.ids.size() == 3);
        assert(flat.vals.size() == 3);

        double sum = 0.0;
        for (double v : flat.vals) sum += v;
        assert(approx(sum, 12.0));

        for (std::size_t i = 0; i < flat.ids.size(); ++i) {
            assert(flat.ids[i].s == obs);
            if (flat.vals[i] == 3.0) {
                assert(approx(flat.ids[i].p.first, 0.0));
                assert(approx(flat.ids[i].p.second, 0.0));
            }
            if (flat.vals[i] == 4.0) {
                assert(approx(flat.ids[i].p.first, 0.1));
                assert(approx(flat.ids[i].p.second, 1.0));
            }
            if (flat.vals[i] == 5.0) {
                assert(approx(flat.ids[i].p.first, 1.0));
                assert(approx(flat.ids[i].p.second, 2.0));
            }
        }
    }

    std::cout << "\nIndexing INTEGRATION tests passed!\n";
    return 0;
}
