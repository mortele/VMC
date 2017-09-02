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
    m_orbital   = new HydrogenOrbital(alpha);
    m_alpha     = alpha;
    m_beta      = beta;
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

    for (int j = 0; j < electron; j++) {
        Electron*       jElectron   = m_system->getElectrons().at(j);
        const double    rij         = computeInterEletronDistance(iElectron, jElectron);
        double          factor      = 1 + m_beta * rij;
        m_jastrowGradient(j, electron) = spinCoefficient(iElectron,jElectron) / (factor*factor);
    }
    for (int j = electron+1; j < m_numberOfElectrons; j++) {
        Electron*       jElectron   = m_system->getElectrons().at(j);
        const double    rij         = computeInterEletronDistance(iElectron, jElectron);
        double          factor      = 1 + m_beta * rij;
        m_jastrowGradient(electron, j) = spinCoefficient(jElectron,iElectron) / (factor*factor);
    }
}


void SlaterWithJastrow::updateJastrowLaplacianTerms(int electron) {
    Electron*   iElectron   = m_system->getElectrons().at(electron);

    for (int j = 0; j < electron; j++) {
        Electron*       jElectron   = m_system->getElectrons().at(j);
        const double    rij         = computeInterEletronDistance(iElectron, jElectron);
        double          factor      = 1 + m_beta * rij;
        m_jastrowLaplacianTerms(j, electron) = -spinCoefficient(iElectron,jElectron) / (factor*factor*factor);
    }
    for (int j = electron+1; j < m_numberOfElectrons; j++) {
        Electron*       jElectron   = m_system->getElectrons().at(j);
        const double    rij         = computeInterEletronDistance(iElectron, jElectron);
        double          factor      = 1 + m_beta * rij;
        m_jastrowLaplacianTerms(electron, j) = -spinCoefficient(jElectron,iElectron) / (factor*factor*factor);
    }
}

void SlaterWithJastrow::computeJastrowLaplacian() {
    double sum = 0;
    for (int k = 0; k < m_numberOfElectrons; k++) {
        Electron* kElectron = m_system->getElectrons().at(k);

        for (int i = 0; i < k; i++) {
            Electron* iElectron = m_system->getElectrons().at(i);
            const double rik = computeInterEletronDistance(kElectron, iElectron);
            sum += 2/rik * m_jastrowGradient(i,k) + m_jastrowLaplacianTerms(i,k);
        }
        for (int i = k+1; i < m_numberOfElectrons; i++) {
            Electron* iElectron = m_system->getElectrons().at(i);
            const double rik = computeInterEletronDistance(iElectron, kElectron);
            sum += 2/rik * m_jastrowGradient(k,i) + m_jastrowLaplacianTerms(k,i);
        }
    }
    m_jastrowLaplacian = -0.5*sum;
}

void SlaterWithJastrow::computeSlaterRatio() {
    double sum = 0;
    Electron* iElectron = m_system->getElectrons().at(m_changedElectron);
    const int    i  = iElectron->getSpinIndex();
    const double xi = iElectron->getPosition().at(0);
    const double yi = iElectron->getPosition().at(1);
    const double zi = iElectron->getPosition().at(2);

    if (m_spinChanged==1) {
        for (int j = 0; j < m_numberOfSpinUpElectrons; j++) {
            sum += (*m_orbital)(xi,yi,zi,j) * m_slaterUp(j,i);
        }
    } else {
        for (int j = 0; j < m_numberOfSpinDownElectrons; j++) {
            sum += (*m_orbital)(xi,yi,zi,j) * m_slaterDown(j,i);
        }
    }
    m_Rsd = sum;
}

void SlaterWithJastrow::computeJastrowRatio() {
    double sum = 0;

    int k = m_changedElectron;
    for (int i = 0; i < k; i++) {
        sum += m_correlationMatrix(i,k) - m_correlationMatrixOld(i,k);
    }
    for (int i = k+1; i < m_numberOfElectrons; i++) {
        sum += m_correlationMatrix(k,i) - m_correlationMatrixOld(k,i);
    }

    m_Rc = exp(sum);
}


void SlaterWithJastrow::computeSlaterLaplacian(int electron) {
    Electron* kElectron = m_system->getElectrons().at(electron);
    const int kSpin     = kElectron->getSpin();
    const int eLimit    = (kSpin==1 ? m_numberOfSpinUpElectrons : m_numberOfSpinDownElectrons);

    double value = 0;
    for (int i = 0; i < eLimit; i++) {
        Electron* iElectron = (kSpin==1 ? m_system->getSpinUpElectrons().at(i) : m_system->getSpinDownElectrons().at(i));
        const double xi = iElectron->getPosition().at(0);
        const double yi = iElectron->getPosition().at(1);
        const double zi = iElectron->getPosition().at(2);

        for (int j = 0; j < eLimit; j++) {
            const double jLaplacian = m_orbital->computeLaplacian(xi,yi,zi,j);
            value += (kSpin==1 ? m_slaterUp(j,i) : m_slaterDown(j,i)) * jLaplacian;
        }
    }
    if (kSpin==1)   m_slaterLaplacianUp = value;
    else            m_slaterLaplacianDown = value;
    m_slaterLaplacian = m_slaterLaplacianUp + m_slaterLaplacianDown;
}

void SlaterWithJastrow::fillCorrelationMatrix() {
    for (int k = 0; k < m_numberOfElectrons; k++) {
        Electron* kElectron = m_system->getElectrons().at(k);

        for (int i = k+1; i < m_numberOfElectrons; i++) {
            Electron* iElectron = m_system->getElectrons().at(i);
            m_correlationMatrix(k,i) = computeJastrowFactor(kElectron,iElectron);
        }
    }
    for (int k = 0; k < m_numberOfElectrons; k++) {
        for (int i = k+1; i < m_numberOfElectrons; i++) {
            m_correlationMatrix(i,k) = m_correlationMatrix(k,i);
        }
    }
}

void SlaterWithJastrow::updateCorrelationsMatrix() {
    int k = m_changedElectron;
    Electron* kElectron = m_system->getElectrons().at(m_changedElectron);
    for (int i = 0; i < m_changedElectron; i++) {
        Electron* iElectron = m_system->getElectrons().at(i);
        m_correlationMatrix(i,k) = computeJastrowFactor(iElectron,kElectron);
        m_correlationMatrix(k,i) = m_correlationMatrix(i,k);
    }
    for (int i = m_changedElectron+1; i < m_numberOfElectrons; i++) {
        Electron* iElectron = m_system->getElectrons().at(i);
        m_correlationMatrix(k,i) = computeJastrowFactor(kElectron,iElectron);
        m_correlationMatrix(i,k) = m_correlationMatrix(k,i);
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
            if (kSpin==1) {
                const int spinIndex = kElectron->getSpinIndex();
                m_quantumForce(electron, dimension) = 2 * m_slaterGradientUp(spinIndex,dimension);
                if (m_jastrow) m_quantumForce(electron,dimension) += 2 * sum;
                m_energyCrossTerm -= 0.5*sum*sum + (m_slaterGradientUp(spinIndex,dimension) * sum);
            } else {
                const int spinIndex = kElectron->getSpinIndex();
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
    for (unsigned long electron = 0; electron < m_system->getElectrons().size(); electron++) {
        m_system->getElectrons().at(electron)->setIndex(electron);
    }
    for (unsigned long electron = 0; electron < m_system->getSpinUpElectrons().size(); electron++) {
        m_system->getSpinUpElectrons().at(electron)->setSpinIndex(electron);
    }
    for (unsigned long electron = 0; electron < m_system->getSpinDownElectrons().size(); electron++) {
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
    m_jastrowGradient       = zeros<mat>(m_numberOfElectrons,m_numberOfElectrons);
    m_jastrowLaplacianTerms = zeros<mat>(m_numberOfElectrons,m_numberOfElectrons);

    for (int electron = 0; electron < m_numberOfElectrons; electron++) {
        updateSlaterGradient(1.0, electron);
        computeSlaterLaplacian(electron);   // This is called N-2 too many times,
                                            // but that shouldnt affect runtime
                                            // in any substantial way.
        if (m_jastrow) updateJastrowGradient(electron);
        if (m_jastrow) updateJastrowLaplacianTerms(electron);
    }

    m_quantumForce = zeros<mat>(m_numberOfElectrons,3);
    computeQuantumForce();

    m_positionsOld = zeros<mat>(m_numberOfElectrons, 3);
    for (int electron = 0; electron < m_numberOfElectrons; electron++) {
        Electron* iElectron = m_system->getElectrons().at(electron);

        for (int dimension = 0; dimension < 3; dimension++) {
            const double x = iElectron->getPosition().at(dimension);
            m_positionsOld(electron,dimension) = x;
        }
    }

    m_correlationMatrix    = zeros<mat>(m_numberOfElectrons, m_numberOfElectrons);
    m_correlationMatrixOld = zeros<mat>(m_numberOfElectrons, m_numberOfElectrons);
    fillCorrelationMatrix();
    m_correlationMatrixOld = m_correlationMatrix;

}

void SlaterWithJastrow::passProposedChangeToWaveFunction(int electronChanged, int dimensionChanged) {
    m_changedElectron   = electronChanged;
    m_changedDimension  = dimensionChanged;
    m_spinChanged       = m_system->getElectrons().at(electronChanged)->getSpin();
    updateCorrelationsMatrix();
}

void SlaterWithJastrow::updateWaveFunctionAfterAcceptedStep() {
    // Update the matrix of old positions.
    Electron* iElectron = m_system->getElectrons().at(m_changedElectron);
    for (int dimension = 0; dimension < 3; dimension++) {
        const double x = iElectron->getPosition().at(dimension);
        m_positionsOld(m_changedElectron, dimension) = x;
    }
    m_correlationMatrixOld = m_correlationMatrix;
}

double SlaterWithJastrow::computeWaveFunctionRatio(int changedElectronIndex) {
    m_changedElectron = changedElectronIndex;
    computeSlaterRatio();
    if (m_jastrow) {
        computeJastrowRatio();
        return m_Rc*m_Rsd;
    } else {
        return m_Rsd;
    }
}

double SlaterWithJastrow::evaluateLaplacian() {
    computeSlaterLaplacian(m_changedElectron);
    computeQuantumForce();
    if (m_jastrow) {
        updateJastrowGradient(m_changedElectron);
        updateJastrowLaplacianTerms(m_changedElectron);
        computeJastrowLaplacian();
    }

    if (m_jastrow) {
        // Multiply the cross-term by (-2) to cancel the -0.5 multiplication in
        // Hamiltonian::computeKineticEnergy(), since this term alone shouldnt
        // have the factor applied to it but the two other terms _should_.
        m_laplacian = m_slaterLaplacian + m_jastrowLaplacian + (-2)*m_energyCrossTerm;
    } else {
        m_laplacian = m_slaterLaplacian;
    }

    return m_laplacian;
}




