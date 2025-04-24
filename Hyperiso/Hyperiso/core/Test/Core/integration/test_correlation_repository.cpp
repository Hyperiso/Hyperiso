#include "CorrelationRepo.h"
#include "Parameter.h"
#include <iostream>
#include <cassert>
#include <cmath>

void test_paramid_correlation_integration() {
    std::cout << "\n-- test_paramid_correlation_integration --" << std::endl;

    // Création d'identifiants de paramètres
    ParamId alphaS("GAUGE", LhaID(3)); // g3
    ParamId MW("MASS", LhaID(24));     // W boson
    ParamId MZ("MASS", LhaID(23));     // Z boson

    // Définir des corrélations entre ces paramètres
    auto matrix = std::make_shared<CorrelationMatrixPair<ParamId>>();
    matrix->emplace({alphaS, MW}, 0.1, 0.05);
    matrix->emplace({alphaS, MZ}, 0.2, 0.1);
    matrix->emplace({MW, MZ}, 0.3, 0.15);

    CorrelationRepository repo;
    repo.set_correlation_matrix(matrix);

    // Vérification de la cohérence des corrélations croisées
    auto c1 = repo.get_correlation(alphaS, MW);
    auto c2 = repo.get_correlation(MZ, alphaS);
    auto c3 = repo.get_correlation(MZ, MW);

    assert(std::abs(c1.first - 0.1) < 1e-6);
    assert(std::abs(c1.second - 0.05) < 1e-6);

    assert(std::abs(c2.first - 0.2) < 1e-6);
    assert(std::abs(c2.second - 0.1) < 1e-6);

    assert(std::abs(c3.first - 0.3) < 1e-6);
    assert(std::abs(c3.second - 0.15) < 1e-6);
}

void test_observables_correlation_integration() {
    std::cout << "\n-- test_observables_correlation_integration --" << std::endl;

    // Définition d'observables physiques
    Observables RD = Observables::R_D;
    Observables RDstar = Observables::R_DSTAR;
    Observables RTau = Observables::R_TAU_NU;

    auto matrix = std::make_shared<CorrelationMatrixPair<Observables>>();
    matrix->emplace({RD, RDstar}, 0.25, 0.05);
    matrix->emplace({RD, RTau}, 0.15, 0.03);

    CorrelationRepository repo;
    repo.set_correlation_matrix(matrix);

    // Vérifications
    auto c1 = repo.get_correlation(RD, RDstar);
    auto c2 = repo.get_correlation(RTau, RD);

    assert(std::abs(c1.first - 0.25) < 1e-6);
    assert(std::abs(c1.second - 0.05) < 1e-6);

    assert(std::abs(c2.first - 0.15) < 1e-6);
    assert(std::abs(c2.second - 0.03) < 1e-6);
}

void test_merge_paramid_and_access_full() {
    std::cout << "\n-- test_merge_paramid_and_access_full --" << std::endl;

    ParamId p1("BLOCK", LhaID(1));
    ParamId p2("BLOCK", LhaID(2));
    ParamId p3("BLOCK", LhaID(3));

    auto m1 = std::make_shared<CorrelationMatrixPair<ParamId>>();
    m1->emplace({p1, p2}, 0.6, 0.4);

    auto m2 = std::make_shared<CorrelationMatrixPair<ParamId>>();
    m2->emplace({p2, p3}, 0.3, 0.2);

    CorrelationRepository repo;
    repo.set_correlation_matrix(m1);
    repo.merge_correlation_matrix(m2);

    auto c1 = repo.get_correlation(p1, p2);
    auto c2 = repo.get_correlation(p2, p3);

    assert(std::abs(c1.first - 0.6) < 1e-6);
    assert(std::abs(c2.second - 0.2) < 1e-6);
}

int main() {
    std::cout << "== Running INTEGRATION tests for CorrelationRepository ==\n";

    test_paramid_correlation_integration();
    test_observables_correlation_integration();
    test_merge_paramid_and_access_full();

    std::cout << "\n✅ All CorrelationRepository integration tests passed!\n" << std::endl;
    return 0;
}
