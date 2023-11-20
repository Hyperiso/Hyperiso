#include "gtest/gtest.h"
#include "RareDecayCalculator.h"

TEST(RareDecayCalculatorTest, CalculatesCorrectBranchingRatio) {
    ModelParameters params(125.0, 0.1); // Exemple de paramètres
    RareDecayCalculator calculator(params);

    double expectedRatio = 12.5; // Valeur attendue pour ces paramètres
    EXPECT_NEAR(calculator.calculateBranchingRatio(), expectedRatio, 0.001);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
