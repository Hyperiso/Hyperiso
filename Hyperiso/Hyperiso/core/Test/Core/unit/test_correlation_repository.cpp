#include "CorrelationRepo.h"
#include "Parameter.h"
#include <iostream>
#include <cassert>

void test_matrix_pair_basic_access() {
    std::cout << "\n-- test_matrix_pair_basic_access --" << std::endl;

    CorrelationMatrixPair<ParamId> matrix;
    
    ParamId p1("BLOCK", LhaID(1));
    ParamId p2("BLOCK", LhaID(2));

    matrix.emplace(std::make_pair(p1, p2), 0.5, 0.2);

    auto corr = matrix.at({p1, p2});
    assert(std::abs(corr.first - 0.5) < 1e-6);
    assert(std::abs(corr.second - 0.2) < 1e-6);

    auto corr_sym = matrix.at({p2, p1});
    assert(std::abs(corr_sym.first - 0.5) < 1e-6);
    assert(std::abs(corr_sym.second - 0.2) < 1e-6);

    auto corr_zero = matrix.at({p1, p1});
    assert(corr_zero.first == 0.0 && corr_zero.second == 0.0);
}

void test_repo_paramid_set_get() {
    std::cout << "\n-- test_repo_paramid_set_get --" << std::endl;

    auto matrix = std::make_shared<CorrelationMatrixPair<ParamId>>();
    ParamId p1("BLOCK", LhaID(1));
    ParamId p2("BLOCK", LhaID(2));
    matrix->emplace({p1, p2}, 0.3, 0.1);

    CorrelationRepository repo;
    repo.set_correlation_matrix(matrix);

    auto corr = repo.get_correlation(p1, p2);
    assert(std::abs(corr.first - 0.3) < 1e-6);
    assert(std::abs(corr.second - 0.1) < 1e-6);

    auto corr_rev = repo.get_correlation(p2, p1);
    assert(std::abs(corr_rev.first - 0.3) < 1e-6);
    assert(std::abs(corr_rev.second - 0.1) < 1e-6);
}

void test_repo_observables_set_get() {
    std::cout << "\n-- test_repo_observables_set_get --" << std::endl;

    auto matrix = std::make_shared<CorrelationMatrixPair<Observables>>();
    matrix->emplace({Observables::R_D, Observables::R_DSTAR}, 0.25, 0.05);

    CorrelationRepository repo;
    repo.set_correlation_matrix(matrix);

    auto corr = repo.get_correlation(Observables::R_D, Observables::R_DSTAR);
    assert(std::abs(corr.first - 0.25) < 1e-6);
    assert(std::abs(corr.second - 0.05) < 1e-6);
}

void test_merge_correlation_matrix() {
    std::cout << "\n-- test_merge_correlation_matrix --" << std::endl;

    CorrelationRepository repo;

    auto m1 = std::make_shared<CorrelationMatrixPair<Observables>>();
    m1->emplace({Observables::R_D, Observables::R_DSTAR}, 0.1, 0.2);
    repo.set_correlation_matrix(m1);

    auto m2 = std::make_shared<CorrelationMatrixPair<Observables>>();
    m2->emplace({Observables::R_TAU_NU, Observables::R_D}, 0.3, 0.4);
    repo.merge_correlation_matrix(m2);

    auto corr1 = repo.get_correlation(Observables::R_D, Observables::R_DSTAR);
    auto corr2 = repo.get_correlation(Observables::R_TAU_NU, Observables::R_D);

    assert(std::abs(corr1.first - 0.1) < 1e-6);
    assert(std::abs(corr2.second - 0.4) < 1e-6);
}

int main() {
    std::cout << "== Running UNIT tests for CorrelationRepository ==\n";

    test_matrix_pair_basic_access();
    test_repo_paramid_set_get();
    test_repo_observables_set_get();
    test_merge_correlation_matrix();

    std::cout << "\n✅ All CorrelationRepository unit tests passed!\n" << std::endl;
    return 0;
}
