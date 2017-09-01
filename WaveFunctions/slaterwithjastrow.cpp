#include "slaterwithjastrow.h"
#include "system.h"
#include "electron.h"
#include <armadillo>
#include <vector>

using arma::mat;
using arma::zeros;
using std::vector;


SlaterWithJastrow::SlaterWithJastrow(System*    system,
                                     double     alpha,
                                     double     beta) :
    WaveFunction(system) {
    m_orbital = new HydrogenOrbital(alpha);
}


void SlaterWithJastrow::updateSlaterGradient(double Rsd, int electron) {
    int spin    = m_system->getElectrons().at(electron)->getSpin();
    int eLimit  = (spin==1 ? m_numberOfSpinUpElectrons : m_numberOfSpinDownElectrons);

    const double x = m_system->getElectrons().at(electron)->getPosition().at(0);
    const double y = m_system->getElectrons().at(electron)->getPosition().at(1);
    const double z = m_system->getElectrons().at(electron)->getPosition().at(2);

    for (int dimension = 0; dimension < 3; dimension++) {
        double sum = 0;
        for (int j = 0; j < eLimit; j++) {
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


void SlaterWithJastrow::updateJastrowLaplacian(int electron) {
    Electron*   iElectron   = m_system->getElectrons().at(electron);
    const int   iSpin       = iElectron->getSpin();

    for (int j = 0; j < electron; j++) {
        Electron*       jElectron   = m_system->getElectrons().at(j);
        const int       jSpin       = jElectron->getSpin();
        const double    rij         = computeInterEletronDistance(iElectron, jElectron);
        double          factor      = 1 + m_beta * rij;
        m_jastrowGradient(j, electron) = -(iSpin==jSpin ? 0.25 : 0.5) / (factor*factor*factor);
    }
    for (int j = electron+1; j < m_numberOfElectrons; j++) {
        Electron*       jElectron   = m_system->getElectrons().at(j);
        const int       jSpin       = jElectron->getSpin();
        const double    rij         = computeInterEletronDistance(iElectron, jElectron);
        double          factor      = 1 + m_beta * rij;
        m_jastrowGradient(electron, j) = -(iSpin==jSpin ? 0.25 : 0.5) / (factor*factor*factor);
    }
}

void SlaterWithJastrow::computeQuantumForce() {
    m_energyCrossTerm = 0;
    for (int electron = 0; electron < m_numberOfElectrons; electron++) {
        Electron*   kElectron   = m_system->getElectrons().at(electron);
        const int   kSpin       = kElectron->getSpin();

        for (int dimension = 0; dimension < 3; dimension++) {
            double sum = 0;
            for (int i = 0; i < electron; i++) {
                Electron*    iElectron  = m_system->getElectrons().at(electron);
                const double xk         = kElectron->getPosition().at(dimension);
                const double xi         = iElectron->getPosition().at(dimension);
                const double rik        = computeInterEletronDistance(kElectron, iElectron);
                sum += (xk - xi) / rik * m_jastrowGradient(i, electron);
            }
            for (int i = electron+1; i < m_numberOfElectrons; i++) {
                Electron*    iElectron  = m_system->getElectrons().at(electron);
                const double xk         = kElectron->getPosition().at(dimension);
                const double xi         = iElectron->getPosition().at(dimension);
                const double rik        = computeInterEletronDistance(iElectron, kElectron);
                sum -= (xi - xk) / rik * m_jastrowGradient(i, electron);
            }
            if (spin==1) {
                const spinIndex = kElectron->getSpinIndex();
                m_quantumForce(electron, dimension) = 2 * m_slaterGradientUp(spinIndex,dimension);
                if (m_jastrow) m_quantumForce(electron,dimension) += 2 * sum;
                m_energyCrossTerm -= 0.5*sum*sum + (m_slaterGradientUp(spinIndex,dimension) * sum);
            } else {
                const spinIndex = kElectron->getSpinIndex();
                m_quantumForce(electron, dimension) = 2 * m_slaterGradientDown(spinIndex,dimension);
                if (m_jastrow) m_quantumForce(electron,dimension) += 2 * sum;
                m_energyCrossTerm -= 0.5*sum*sum + (m_slaterGradientDown(spinIndex,dimension) * sum);
            }
        }
    }
}


void SlaterWithJastrow::evaluateWaveFunctionInitial() {
    const int eUp   = m_numberOfSpinUpElectrons;
    const int eDown = m_numberOfSpinDownElectrons;

    // Label all the electrons so we know which index in the total electron
    // array correspond to the same electrons in the spin-restricted electron
    // arrays.
    for (int electron = 0; electron < m_system->getElectrons().size(); electron++) {
        m_system->getElectrons().at(electron)->setIndex(electron);
    }
    for (int electron = 0; electron < m_system->getSpinUpElectrons().size(); electron++) {
        m_system->getSpinUpElectrons().at(electron)->setSpinIndex(electron);
    }
    for (int electron = 0; electron < m_system->getSpinDownElectrons().size(); electron++) {
        m_system->getSpinDownElectrons().at(electron)->setSpinIndex(electron);
    }


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
    m_jastrowGradient    = zeros<mat>(m_numberOfElectrons,m_numberOfElectrons);

    for (int electron = 0; electron < m_numberOfElectrons; electron++) {
        updateSlaterGradient(1.0, electron);
        if (m_jastrow) updateJastrowGradient(electron);
        if (m_jastrow) updateJastrowLaplacian(electron);
    }

    m_quantumForce = zeros<mat>(m_numberOfElectrons,3);
}

void SlaterWithJastrow::updateWaveFunction(int electronChanged, int dimensionChanged) {
    m_changedElectron   = electronChanged;
    m_changedDimension  = dimensionChanged;
    m_spinChanged       = m_system->getElectrons().at(electronChanged)->getSpin();
}




