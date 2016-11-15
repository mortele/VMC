#include "hartreefockbasisparser.h"
#include "WaveFunctions/Orbitals/primitivegaussian.h"
#include "WaveFunctions/Orbitals/contractedgaussian.h"
#include <fstream>

using std::cout;
using std::endl;
using std::string;
using arma::vec;
using arma::mat;
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
 * cUp00    cUp01   cUp02   cUp03   cUp04   ...
 * cUp10    cUp11   cUp12   cUp13   cUp14   ...
 * cUp20    cUp21   cUp22   cUp23   cUp24   ...
 * ...
 * cDown00  cDown01 cDown02 cDown03 cDown04 ...
 * cDown10  cDown11 cDown12 cDown13 cDown14 ...
 * cDown20  cDown21 cDown22 cDown23 cDown24 ...
 * ...
 * ____________________________________________________________________________
 */
void HartreeFockBasisParser::parseBasisFile(string fileName) {
    std::ifstream basisFile;
    basisFile.open(fileName, std::ios::in);
    if (! basisFile.is_open()) {
        std::cout << "Unable to open file: " << fileName << std::endl;
        return;
    }

    basisFile >> m_numberOfElectrons
              >> m_numberOfSpinUpElectrons
              >> m_numberOfSpinDownElectrons
              >> m_basisSize
              >> m_numberOfAtoms;

    m_atomCharges  .reserve(m_numberOfAtoms);
    m_atomPositions.reserve(m_numberOfAtoms);

    for (int atom = 0; atom < m_numberOfAtoms; atom++) {
        int    Z = 0;
        double x = 0;
        double y = 0;
        double z = 0;
        basisFile >> Z >> x >> y >> z;
        m_atomCharges.push_back(Z);
        m_atomPositions.push_back(vec{x,y,z});
    }
    m_basis.reserve(m_basisSize);

    for (int basisFunction = 0; basisFunction < m_basisSize; basisFunction++) {
        int primitives;
        double atomx, atomy, atomz;
        basisFile >> primitives >> atomx >> atomy >> atomz;

        m_basis.push_back(new ContractedGaussian(primitives, atomx, atomy, atomz));

        for (int primitive = 0; primitive < primitives; primitive++) {
            int    x = 0;
            int    y = 0;
            int    z = 0;
            double exponent = 0;
            double constant = 0;
            basisFile >> x >> y >> z >> exponent >> constant;

            m_basis.at(basisFunction)->addPrimitive(new PrimitiveGaussian(x, y, z, exponent, constant));
        }
    }

    m_spinUpCoefficients   = zeros<mat>(m_basisSize, m_numberOfSpinUpElectrons);
    m_spinDownCoefficients = zeros<mat>(m_basisSize, m_numberOfSpinDownElectrons);
    for (int electron = 0; electron < m_numberOfSpinUpElectrons; electron++) {
        for (int basisFunction = 0; basisFunction < m_basisSize; basisFunction++) {
            double coefficient = 0;
            basisFile >> coefficient;
            m_spinUpCoefficients(basisFunction, electron) = coefficient;
        }
    }
    for (int electron = 0; electron < m_numberOfSpinDownElectrons; electron++) {
        for (int basisFunction = 0; basisFunction < m_basisSize; basisFunction++) {
            double coefficient = 0;
            basisFile >> coefficient;
            m_spinUpCoefficients(basisFunction, electron) = coefficient;
        }
    }
    basisFile.close();
}





