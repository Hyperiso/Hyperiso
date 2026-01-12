#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <unordered_map>

#include "Block.h"
#include "Parameter.h"
#include "SourcesView.h"
#include "Include.h"

static std::shared_ptr<Parameter> mkp(const std::string& blk, int code, double v) {
    return std::make_shared<Parameter>(ParamId{ParameterType::SM, blk, code}, v, 0.0, 0.0);
}

int main() {
    std::cout << "== Running Block integration tests ==\n";

    LhaID k(100);

    auto src1 = std::make_shared<Block>(); src1->blockname = "SRC1"; src1->bind_self(src1);
    auto src2 = std::make_shared<Block>(); src2->blockname = "SRC2"; src2->bind_self(src2);

    src1->store(k, mkp("SRC1", 100, 1.0));
    src2->store(k, mkp("SRC2", 100, 2.0));

    int runsSum = 0;
    auto depSum = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{ {"SRC1", src1}, {"SRC2", src2} },
        [k, &runsSum](const auto& blocks, std::shared_ptr<DependentBlock> self) {
            ++runsSum;
            double v = blocks.get_val("SRC1", k) + blocks.get_val("SRC2", k);
            if (!self->contains(k)) self->store(k, mkp("SUM", 100, 0.0));
            self->assign(k, v);
        }
    );
    depSum->blockname = "SUM";
    depSum->bind_self(depSum);
    depSum->init();

    int runsPost = 0;
    auto depPost = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{ {"SUM", depSum} },
        [k, &runsPost](const auto& blocks, std::shared_ptr<DependentBlock> self) {
            ++runsPost;
            double base = blocks.get_val("SUM", k);
            if (!self->contains(k)) self->store(k, mkp("POST", 100, 0.0));
            self->assign(k, 3.0 * base);
        }
    );
    depPost->blockname = "POST";
    depPost->bind_self(depPost);
    depPost->init();

    int runsLeaf = 0;
    auto leaf = std::make_shared<DependentBlock>(
        std::unordered_map<std::string, std::shared_ptr<Block>>{ {"POST", depPost} },
        [k, &runsLeaf](const auto& blocks, std::shared_ptr<DependentBlock> self) {
            ++runsLeaf;
            double x = blocks.get_val("POST", k);
            if (!self->contains(k)) self->store(k, mkp("LEAF", 100, 0.0));
            self->assign(k, x + 1.0);
        }
    );
    leaf->blockname = "LEAF";
    leaf->bind_self(leaf);
    leaf->init();

    assert(runsSum == 0 && runsPost == 0 && runsLeaf == 0);

    assert(std::abs(leaf->retrieve(k)->get_val() - (3.0 * (1.0 + 2.0) + 1.0)) < 1e-12);
    assert(runsSum == 1 && runsPost == 1 && runsLeaf == 1);

    src2->assign(k, 5.0);
    assert(runsSum == 1 && runsPost == 1 && runsLeaf == 1);

    assert(std::abs(leaf->retrieve(k)->get_val() - (3.0 * (1.0 + 5.0) + 1.0)) < 1e-12);
    assert(runsSum == 2 && runsPost == 2 && runsLeaf == 2);

    depSum->freeze();
    src1->assign(k, 10.0);

    assert(std::abs(leaf->retrieve(k)->get_val() - (3.0 * (1.0 + 5.0) + 1.0)) < 1e-12);

    depSum->unfreeze();
    assert(std::abs(leaf->retrieve(k)->get_val() - (3.0 * (10.0 + 5.0) + 1.0)) < 1e-12);

    depPost->clear_above();
    src2->assign(k, 7.0);

    assert(std::abs(depSum->retrieve(k)->get_val() - (10.0 + 7.0)) < 1e-12);

    assert(std::abs(depPost->retrieve(k)->get_val() - 45.0) < 1e-12);
    assert(std::abs(leaf->retrieve(k)->get_val() - 46.0) < 1e-12);

    depPost->addObserver(leaf);
    depPost->clear_below();

    src1->assign(k, 0.0);

    assert(std::abs(depSum->retrieve(k)->get_val() - (0.0 + 7.0)) < 1e-12); 
    assert(std::abs(depPost->retrieve(k)->get_val() - 45.0) < 1e-12);
    assert(std::abs(leaf->retrieve(k)->get_val() - 46.0) < 1e-12);

    std::cout << "\n Block integration suite passed.\n";
    return 0;
}
