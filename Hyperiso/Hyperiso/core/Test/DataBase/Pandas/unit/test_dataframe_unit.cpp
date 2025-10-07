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
    std::cout << "== DataFrame UNIT ==\n";

    DataFrame df;

    df.addColumn<int>("A");
    df.addColumn<double>("B");
    df.addColumn<std::string>("C");

    for (int i = 0; i < 4; ++i) {
        df.addValueToColumn<int>("A", i+1);
        df.addValueToColumn<double>("B", 0.5 * (i+1));
        df.addValueToColumn<std::string>("C", std::string("s") + std::to_string(i+1));
    }

    assert(df.shape[0] == 4);
    assert(df.shape[1] == 3);
    assert(df.getRowCount() == 4);

    assert(df.iat<int>(0, "A") == 1);
    assert(std::abs(df.iat<double>(3, "B") - 2.0) < 1e-12);
    assert(df.iat<std::string>(1, "C") == "s2");

    df.setIndex({"r1","r2","r3","r4"});
    assert(df.at<int>("r3","A") == 3);
    assert(df.at<std::string>("r4","C") == "s4");

    const auto& idx = df.getIndex();
    assert(idx.size() == 4 && idx[0] == "r1");

    df.operator[]<double>("D").add(10.0);
    df.addValueToColumn<double>("D", 11.0);
    df.addValueToColumn<double>("D", 12.0);
    df.addValueToColumn<double>("D", 13.0);
    assert(df.shape[1] == 4);
    assert(df.getRowCount() == 4);
    assert(std::abs(df.iat<double>(0, "D") - 10.0) < 1e-12);

    auto h2 = df.head(2);
    assert(h2.getRowCount() == 2);
    assert(h2.iat<int>(1, "A") == 2);
    auto t3 = df.tail(3);
    assert(t3.getRowCount() == 3);
    assert(t3.iat<std::string>(2, "C") == "s4");

    df.describe();

    {
        CSVOptions opt;
        opt.hasIndex = false;
        df._set_csv_options(opt);

        auto tmp = fs::temp_directory_path() / "df_unit_no_index.csv";
        df.to_csv(tmp.string());
        auto content = slurp(tmp);

        assert(content.find("A,B,C,D\n") == 0);
        assert(content.find("1,0.5,s1,10") != std::string::npos);
    }

    {
        CSVOptions opt;
        opt.hasIndex = true;
        df._set_csv_options(opt);

        auto tmp = fs::temp_directory_path() / "df_unit_with_index.csv";
        df.to_csv(tmp.string());
        auto content = slurp(tmp);

        assert(content.find("Index,A,B,C,D\n") == 0);
        assert(content.find("0,1,0.5,s1,10") != std::string::npos);
        assert(content.find(",\n") == std::string::npos || content.find(",\n") > 10);
    }

    {
        DataFrame df2;
        df2.addColumn<int>("X");
        df2.addValueToColumn<int>("X", 42);
        bool threw = false;
        try { (void)df2.getIndex(); } catch (const std::runtime_error&) { threw = true; }
        assert(threw);
    }

    std::cout << "✅ UNIT OK\n";
    return 0;
}
