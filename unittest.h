#pragma once

class UnitTest {
private:
    static class System* setupNewTestSystem();

public:
    static bool runAllTests();
    static bool testHydrogen();
    static bool testHelium();
    static bool testHeliumWithJastrow();
    static bool testNonInteractingHelium();
    static bool testNumericalLaplacian();
};

