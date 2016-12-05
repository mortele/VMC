#pragma once

class UnitTest {
private:
    static class System* setupNewTestSystem();

public:
    static bool runAllTests();
    static bool testHydrogen();
    static bool testHelium();
    static bool testHeliumWithJastrowNumerical();
    static bool testNonInteractingHelium();
    static bool testNumericalLaplacian();
    static bool testDirectSlaterHelium();
    static bool testDirectSlaterBeryllium();
    static bool testDirectSlaterWithJastrowHelium();
    static bool testDirectSlaterWithJastrowBeryllium();
    static bool testGaussianSlaterHydrogenMolecule();
    static bool HO3d();
};

