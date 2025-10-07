#include <cassert>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include "DataFrame.h"
#include "CSVOptions.h"

namespace fs = std::filesystem;

static std::string slurp(const fs::path& p) {
    std::ifstream f(p);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

int main() {
    std::cout << "== DataFrame INTEGRATION ==\n";

    DataFrame df;

    df.addColumn<int>("id");
    df.addColumn<double>("x");
    df.addColumn<std::string>("label");

    for (int i = 0; i < 6; ++i) {
        df.addValueToColumn<int>("id", 100 + i);
        df.addValueToColumn<double>("x", 1.25 * i);
        df.addValueToColumn<std::string>("label", std::string("L") + std::to_string(i));
    }
    df.setIndex({"i0","i1","i2","i3","i4","i5"});

    assert(df.getRowCount() == 6);
    assert(df.shape[1] == 3);
    assert(df.at<int>("i3", "id") == 103);
    assert(std::abs(df.at<double>("i5", "x") - 6.25) < 1e-12);

    auto h = df.head(3);
    assert(h.getRowCount() == 3);
    assert(h.iat<std::string>(2,"label") == "L2");

    auto t = df.tail(2);
    assert(t.getRowCount() == 2);
    assert(t.iat<int>(0,"id") == 104);
    assert(t.iat<std::string>(1,"label") == "L5");

    CSVOptions opt; opt.hasIndex = true;
    df._set_csv_options(opt);
    auto tmp = fs::temp_directory_path() / "df_integ.csv";
    df.to_csv(tmp.string());
    auto content = slurp(tmp);

    assert(content.find("Index,id,x,label\n") == 0);

    assert(content.find("0,100,0,L0\n") != std::string::npos);
    assert(content.find("5,105,6.25,L5\n") != std::string::npos);

    auto xcol = df.getColumn<double>("x");
    assert(xcol.size() == 6);
    assert(std::abs(xcol.iat(4) - 5.0) < 1e-12);

    std::cout << "✅ INTEGRATION OK\n";
    return 0;
}
