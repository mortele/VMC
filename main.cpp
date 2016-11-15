#include "unittest.h"
#include "hartreefockbasisparser.h"

int main(int, char**) {
    //UnitTest::runAllTests();

    HartreeFockBasisParser parser;
    parser.parseBasisFile("../../HartreeFock/HartreeFockBases/basis-2016-11-15-10.02.27");
}

