#include "slaterwithjastrow.h"
#include "system.h"
#include "electron.h"
#include <armadillo>
#include <vector>

using arma::mat;
using arma::zeros;
using std::vector;


void SlaterWithJastrow::updateSlaterGradient(double Rsd, int electron) {
    int spin = m_system->getElectrons().at(electron)->getSpin();

    const double x = m_system->getElectrons.at(electron)->getPosition().at(0);
    const double y = m_system->getElectrons.at(electron)->getPosition().at(1);
    const double z = m_system->getElectrons.at(electron)->getPosition().at(2);

    for (int dimension = 0; dimension < 3; dimension++) {
        double sum = 0;
        for (int j = 0; j < eUp; j++) {
            sum += m_orbital->computeDerivative(x,y,z,j,dimension);
        }
        if (spin==1) {
            m_slaterGradientUp(electron, dimension) = sum / Rsd;
        } else {
            m_slaterGradientDown(electron, dimension) = sum / Rsd;
        }
    }
}

void SlaterWithJastrow::updateJastrowGradient(int electron) {
    Electron*   iElectron   = m_system->getElectrons().at(electron);
    const int   iSpin       = iElectron->getSpin();

    for (int j = 0; j < electron; j++) {
        Electron*       jElectron   = m_system->getElectrons().at(j);
        const int       jSpin       = jElectron->getSpin();
        const double    rij         = computeInterEletronDistance(iElectron, jElectron);
        double          factor      = 1 + m_beta * rij;
        m_jastrowGradient(j, electron) = (iSpin==jSpin ? 0.25 : 0.5) / (factor*factor);
    }
    for (int j = electron+1; j < m_numberOfElectrons; j++) {
        Electron*       jElectron   = m_system->getElectrons().at(j);
        const int       jSpin       = jElectron->getSpin();
        const double    rij         = computeInterEletronDistance(iElectron, jElectron);
        double          factor      = 1 + m_beta * rij;
        m_jastrowGradient(electron, j) = (iSpin==jSpin ? 0.25 : 0.5) / (factor*factor);
    }
}

SlaterWithJastrow::SlaterWithJastrow(System*    system,
                                     double     alpha,
                                     double     beta) :
    WaveFunction(system) {
    m_orbital = new HydrogenOrbital(alpha);
}

void SlaterWithJastrow::evaluateWaveFunctionInitial() {
    const int eUp   = m_numberOfSpinUpElectrons;
    const int eDown = m_numberOfSpinDownElectrons;

    m_slaterUp   = zeros<mat>(eUp,   eUp);
    m_slaterDown = zeros<mat>(eDown, eDown);

    vector<Electron*> spinUpElectrons   = m_system->getSpinUpElectrons();
    vector<Electron*> spinDownElectrons = m_system->getSpinDownElectrons();

    for (int electron = 0; electron < eUp; electron++) {
        const double x = spinUpElectrons.at(electron)->getPosition().at(0);
        const double y = spinUpElectrons.at(electron)->getPosition().at(1);
        const double z = spinUpElectrons.at(electron)->getPosition().at(2);

        for (int basis = 0; basis < eUp; basis++) {
            m_slaterUp(electron, basis) = (*m_orbital)(x,y,z,basis);
        }
    }
    for (int electron = 0; electron < eDown; electron++) {
        const double x = spinDownElectrons.at(electron)->getPosition().at(0);
        const double y = spinDownElectrons.at(electron)->getPosition().at(1);
        const double z = spinDownElectrons.at(electron)->getPosition().at(2);

        for (int basis = 0; basis < eUp; basis++) {
            m_slaterDown(electron, basis) = (*m_orbital)(x,y,z,basis);
        }
    }
    m_slaterUp      = m_slaterUp.i();
    m_slaterDown    = m_slaterDown.i();

    m_slaterGradientUp   = zeros<mat>(eUp,   3);
    m_slaterGradientDown = zeros<mat>(eDown, 3);

    for (int electron = 0; electron < m_numberOfElectrons; electron++)
        updateSlaterGradient(1.0, electron);
        updateJastrowGradient
    }
}

void SlaterWithJastrow::updateWaveFunction(int electronChanged, int dimensionChanged) {
    m_changedElectron   = electronChanged;
    m_changedDimension  = dimensionChanged;
    m_spinChanged       = m_system->getElectrons().at(electronChanged)->getSpin();
}




