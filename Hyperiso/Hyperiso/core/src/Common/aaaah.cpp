#include <iostream>
#include <cassert>
#include <string>
#include <unordered_set>
#include <set>
#include "General.h" // Replace with your actual header path

int main() {
    LOG_INFO("Testing default constructor...");
    BlockName b1;
    assert(b1.get_alias().empty());

    LOG_INFO("Testing constructor with const char*...");
    BlockName b2("MainBlock");
    assert(b2.hasAlias("MainBlock"));

    LOG_INFO("Testing constructor with std::string...");
    BlockName b3(std::string("BlockA"));
    assert(b3.hasAlias("BlockA"));

    LOG_INFO("Testing constructor with initializer_list...");
    BlockName b4({"Alias1", "Alias2"});
    assert(b4.hasAlias("Alias1") && b4.hasAlias("Alias2"));

    LOG_INFO("Testing constructor with unordered_set...");
    std::unordered_set<std::string> set_aliases = {"X", "Y"};
    BlockName b5(set_aliases);
    assert(b5.hasAlias("X") && b5.hasAlias("Y"));

    LOG_INFO("Testing addAlias...");
    b1.addAlias("Added");
    assert(b1.hasAlias("Added"));

    LOG_INFO("Testing casting to string (should log warning if multiple aliases)...");
    std::string str_from_b4 = static_cast<std::string>(b4);
    assert(!str_from_b4.empty());

    LOG_INFO("Testing to_string()...");
    std::string str2 = b4.to_string();
    assert(!str2.empty());

    LOG_INFO("Testing BlockName == BlockName...");
    BlockName b6({"Alias2", "SomethingElse"});
    assert(b4 == b6);

    LOG_INFO("Testing BlockName != BlockName...");
    assert(!(b4 != b6));

    LOG_INFO("Testing BlockName == std::string...");
    assert(b4 == "Alias1");

    LOG_INFO("Testing std::string == BlockName...");
    assert("Alias2" == b4);

    LOG_INFO("Testing BlockName != std::string...");
    assert(b4 != "NotThere");

    LOG_INFO("Testing std::string + BlockName...");
    BlockName b7 = std::string("pre_") + b4;
    assert(b7.hasAlias("pre_Alias1") && b7.hasAlias("pre_Alias2"));

    LOG_INFO("Testing operator<<...");
    std::cout << "b4: " << b4 << std::endl;

    LOG_INFO("Testing to_upper()...");
    b4.to_upper();
    assert(b4.hasAlias("ALIAS1") && b4.hasAlias("ALIAS2"));

    LOG_INFO("Testing operator<...");
    BlockName b8({"alpha", "beta"});
    BlockName b9({"gamma"});
    assert((b8 < b9) || (b9 < b8) || (b8.get_alias() == b9.get_alias()));

    LOG_INFO("Testing unordered_set with BlockName and hash specialization...");
    std::unordered_set<BlockName> name_set;
    name_set.insert(b2);
    name_set.insert(b3);
    name_set.insert(b4);
    assert(name_set.find(b2) != name_set.end());
    assert(name_set.find(b3) != name_set.end());

    LOG_INFO("All tests passed.");
    return 0;
}
