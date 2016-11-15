#include "hartreefockbasisparser.h"
#include "WaveFunctions/wavefunction.h"
#include <fstream>


using std::cout;
using std::endl;
using std::string;
using arma::vec;
using arma::zeros;


/*
 * File format:
 * ____________________________________________________________________________
 * nElectrons nSpinUp nSpinDown nBasis nAtoms
 * atom0_x atom0_y atom0_z
 * atom1_x atom1_y atom1_z
 * atom2_x atom2_y atom2_z
 * ...
 * nPrimitives_basis0
 * center0_x center0_y center0_z
 * prim0_xExp prim0_yExp prim0_zExp prim0_exp prim0_const
 * prim1_xExp prim1_yExp prim1_zExp prim1_exp prim1_const
 * ...
 * nPrimitives_basis1
 * center1_x center1_y center1_z
 * prim0_xExp prim0_yExp prim0_zExp prim0_exp prim0_const
 * prim1_xExp prim1_yExp prim1_zExp prim1_exp prim1_const
 * ...
 * ____________________________________________________________________________
 */
WaveFunction* HartreeFockBasisParser::parseBasisFile(string fileName) {
    std::ifstream basisFile;
    basisFile.open(fileName, std::ios::in);

    basisFile >> m_numberOfElectrons
              >> m_numberOfSpinUpElectrons
              >> m_numberOfSpinDownElectrons
              >> m_basisSize
              >> m_numberOfAtoms;

    cout << m_numberOfElectrons << " "
         << m_numberOfSpinUpElectrons << " "
         << m_numberOfSpinDownElectrons << " "
         << m_basisSize << " "
         << m_numberOfAtoms << " " << endl;

    for (int atom = 0; atom < m_numberOfAtoms; atom++) {
        double x,y,z;
        basisFile >> x >> y >> z;
        m_atomPositions.push_back(vec{x,y,z});
        cout << x << " " << y << " " << z << endl;
    }

    for (int basisFunction = 0; basisFunction < m_basisSize; basisFunction++) {
        int primitives;
        double atomx, atomy, atomz;
        basisFile >> primitives
                  >> atomx
                  >> atomy
                  >> atomz;

        for (int primitive = 0; primitive < primitives; primitive++) {
            int    x,y,z;
            double exponent, constant;
            basisFile >> x >> y >> z >> exponent >> constant;
            cout << x << " " << y << " " << z << " " << exponent << " " << constant << endl;
        }
    }


    basisFile.close();
}
