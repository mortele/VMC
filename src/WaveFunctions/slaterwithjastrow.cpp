#include "slaterwithjastrow.h"
#include "system.h"
#include "electron.h"
#include "hamiltonian.h"
#include "metropolis.h"
#include <armadillo>
#include <vector>
#include <iomanip>

using arma::mat;
using arma::vec;
using arma::zeros;
using std::vector;
using std::cout;
using std::endl;
using std::setprecision;


SlaterWithJastrow::SlaterWithJastrow(System*    system,
                                     double     beta,
                                     bool       useJastrow) :
    WaveFunction(system) {
    m_beta                      = beta;
    m_useNumericalDerivatives   = false;
    m_jastrow                   = useJastrow;
    m_containsJastrow           = true;
}


void SlaterWithJastrow::updateSlaterGradient(double Rsd, int electron) {
    Electron*   iElectron   = m_system->getElectrons().at(electron);
    int cs              = iElectron->getSpin();
    int i               = iElectron->getSpinIndex();
    int eLimit          = (cs == 1 ? m_numberOfSpinUpElectrons :
                                     m_numberOfSpinDownElectrons);
    mat& slaterInverse  = (cs == 1 ? m_slaterUp         : m_slaterDown);
    mat& slaterGradient = (cs == 1 ? m_slaterGradientUp : m_slaterGradientDown);

    const double x = iElectron->getPosition().at(0);
    const double y = iElectron->getPosition().at(1);
    const double z = iElectron->getPosition().at(2);

    for (int dimension = 0; dimension < m_numberOfDimensions; dimension++) {
        double sum = 0;
        for (int j = 0; j < eLimit; j++) {
            sum += m_orbital->computeDerivative(x,y,z,j,dimension,cs) *
                              slaterInverse(j,i);
        }
        slaterGradient(i, dimension) = sum / Rsd;
    }
}

void SlaterWithJastrow::updateJastrowGradient(int k) {
    mat& a = m_spinMatrix;
    mat& R = m_interElectronDistances;

    for (int i = 0; i < k; i++) {
        const double    rik         = R(k,i);
        double          factor      = 1 + m_beta * rik;
        m_jastrowGradient(i, k) = a(k,i) / (factor*factor);
    }
    for (int j = k+1; j < m_numberOfElectrons; j++) {
        const double    rkj         = R(k,j);
        double          factor      = 1 + m_beta * rkj;
        m_jastrowGradient(k, j) = a(j,k) / (factor*factor);
    }
}


void SlaterWithJastrow::updateJastrowLaplacianTerms(int k) {
    mat& a          = m_spinMatrix;
    mat& R          = m_interElectronDistances;
    mat& laplacianJ = m_jastrowLaplacianTerms;

    for (int j = 0; j < k; j++) {
        const double    rkj           = R(k,j);
        double          factor        = 1 + m_beta * rkj;
        laplacianJ(j, k) = -2*a(j,k) *m_beta / (factor*factor*factor);
    }
    for (int j = k+1; j < m_numberOfElectrons; j++) {
        const double    rkj           = R(k,j);
        double          factor        = 1 + m_beta * rkj;
        laplacianJ(k, j) = -2*a(k,j)*m_beta / (factor*factor*factor);
    }
}

void SlaterWithJastrow::computeJastrowLaplacian() {
    mat& R          = m_interElectronDistances;
    mat& laplacianJ = m_jastrowLaplacianTerms;
    mat& gradientJ  = m_jastrowGradient;

    double sum = 0;
    for (int k = 0; k < m_numberOfElectrons; k++) {
        for (int i = 0; i < k; i++) {
            const double rik = R(i,k);
            sum += 2/rik * gradientJ(i,k) + laplacianJ(i,k);
        }
        for (int i = k+1; i < m_numberOfElectrons; i++) {
            const double rik = R(i,k);
            sum += 2/rik * gradientJ(k,i) + laplacianJ(k,i);
        }
    }
    m_jastrowLaplacian = sum;
}

void SlaterWithJastrow::updateSlaterInverse() {
    Electron* iElectron = m_system->getElectrons().at(m_changedElectron);
    int iSpin   = iElectron->getSpin();
    int i       = iElectron->getSpinIndex();
    const double x = iElectron->getPosition().at(0);
    const double y = iElectron->getPosition().at(1);
    const double z = iElectron->getPosition().at(2);

    int  eLimit;
    mat& newS   = (iSpin==1 ? m_slaterUp    : m_slaterDown);
    mat& oldS   = (iSpin==1 ? m_slaterUpOld : m_slaterDownOld);
    eLimit      = (iSpin==1 ? m_numberOfSpinUpElectrons : m_numberOfSpinDownElectrons);

    for (int k = 0; k < eLimit; k++) {
        for (int j = 0; j < eLimit; j++) {
            if (j != i) {
                double sum = 0;
                for (int l = 0; l < eLimit; l++) {
                    // Maybe should be using the spinIndex instead of global
                    // index for the x,y,z here
                    sum += oldS(l,j) * m_orbital->evaluate(x,y,z,l,iSpin);
                }
                newS(k,j) = oldS(k,j) - oldS(k,i) * sum / m_Rsd;
            } else {
                newS(k,j) = oldS(k,i) / m_Rsd;
            }
        }
    }
}

void SlaterWithJastrow::computeSlaterRatio() {
    Electron* iElectron = m_system->getElectrons().at(m_changedElectron);
    int    i  = iElectron->getSpinIndex();
    double xi = iElectron->getPosition().at(0);
    double yi = iElectron->getPosition().at(1);
    double zi = iElectron->getPosition().at(2);
    double  nElectrons  = (m_spinChanged == 1 ? m_numberOfSpinUpElectrons :
                                                m_numberOfSpinDownElectrons);
    mat&    slater      = (m_spinChanged == 1 ? m_slaterUp : m_slaterDown);


    double sum = 0;
    for (int j = 0; j < nElectrons; j++) {
        sum += m_orbital->evaluate(xi,yi,zi,j,m_spinChanged) * slater(j,i);
    }
    m_Rsd = sum;
    /*if (m_spinChanged==1) {
        for (int j = 0; j < m_numberOfSpinUpElectrons; j++) {
            sum += m_orbital->evaluate(xi,yi,zi,j,1) * m_slaterUp(j,i);
        }
    } else {
        for (int j = 0; j < m_numberOfSpinDownElectrons; j++) {
            sum += m_orbital->evaluate(xi,yi,zi,j,0) * m_slaterDown(j,i);
        }
    }*/
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
    Electron* kElectron  = m_system->getElectrons().at(electron);
    int kSpin      = kElectron->getSpin();
    int nElectrons = (kSpin == 1 ? m_numberOfSpinUpElectrons :
                                   m_numberOfSpinDownElectrons);
    const vector<Electron*>& electrons = (kSpin==1 ?
                                              m_system->getSpinUpElectrons() :
                                              m_system->getSpinDownElectrons());
    const mat& slater = (kSpin==1 ? m_slaterUp : m_slaterDown);
    double& slaterLaplacian = (kSpin==1 ? m_slaterLaplacianUp :
                                          m_slaterLaplacianDown);
    double value = 0;
    for (int i = 0; i < nElectrons; i++) {
        Electron* iElectron = electrons.at(i);
        const double xi = iElectron->getPosition().at(0);
        const double yi = iElectron->getPosition().at(1);
        const double zi = iElectron->getPosition().at(2);

        for (int j = 0; j < nElectrons; j++) {
            double jLaplacian = m_orbital->computeLaplacian(xi,yi,zi,j,kSpin);
            value +=  slater(j,i) * jLaplacian;
        }
    }
    slaterLaplacian = value;
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
    mat& R = m_interElectronDistances;

    int k = m_changedElectron;
    //Electron* kElectron = m_system->getElectrons().at(k);
    for (int i = 0; i < k; i++) {
        //Electron* iElectron = m_system->getElectrons().at(i);
        m_correlationMatrix(i,k) = m_spinMatrix(i,k) * R(i,k) / (1 + m_beta * R(i,k)); //computeJastrowFactor(iElectron,kElectron);
        m_correlationMatrix(k,i) = m_correlationMatrix(i,k);
    }
    for (int i = k+1; i < m_numberOfElectrons; i++) {
        //Electron* iElectron = m_system->getElectrons().at(i);
        m_correlationMatrix(k,i) = m_spinMatrix(k,i) * R(k,i) / (1 + m_beta * R(k,i)); //computeJastrowFactor(kElectron,iElectron);
        m_correlationMatrix(i,k) = m_correlationMatrix(k,i);
    }
}

void SlaterWithJastrow::updateElectronDistanceMatrices() {
    int k = m_changedElectron;
    Electron* kElectron = m_system->getElectrons().at(k);
    mat& r = m_electronPositions;
    mat& R = m_interElectronDistances;

    for (int dimension = 0; dimension < 3; dimension++) {
        r(k,dimension) = kElectron->getPosition().at(dimension);
    }
    for (int i = 0; i < k; i++) {
        double x = r(k,0) - r(i,0);
        double y = r(k,1) - r(i,1);
        double z = r(k,2) - r(i,2);
        R(k,i) = sqrt(x*x + y*y + z*z);
        R(i,k) = R(k,i);
    }
    for (int i = k+1; i < m_numberOfElectrons; i++) {
        double x = r(k,0) - r(i,0);
        double y = r(k,1) - r(i,1);
        double z = r(k,2) - r(i,2);
        R(k,i) = sqrt(x*x + y*y + z*z);
        R(i,k) = R(k,i);
    }
    double x = r(k,0);
    double y = r(k,1);
    double z = r(k,2);
    R(k,k) = sqrt(x*x + y*y + z*z);
}

double SlaterWithJastrow::computeBetaDerivative() {
    mat& a = m_spinMatrix;
    mat& R = m_interElectronDistances;

    double m_betaDerivative = 0;
    for (int i = 0; i < m_numberOfElectrons; i++) {
        for (int j = i+1; j < m_numberOfElectrons; j++) {
            double factor = 1.0 + m_beta * R(i,j);
            m_betaDerivative -= a(i,j) * R(i,j) * R(i,j) / (factor*factor);
        }
    }
    return m_betaDerivative;
}

void SlaterWithJastrow::computeQuantumForce() {
    mat& R = m_interElectronDistances;
    mat& r = m_electronPositions;

    m_energyCrossTerm = 0;
    for (int k = 0; k < m_numberOfElectrons; k++) {
        const int   kSpin       = m_system->getElectrons().at(k)->getSpin();
        const int   kSpinIndex  = m_system->getElectrons().at(k)->getSpinIndex();

        mat& slaterGradient = (kSpin==1 ? m_slaterGradientUp : m_slaterGradientDown);

        for (int j = 0; j < 3; j++) {
            double sum = 0;
            for (int i = 0; i < k ; i++) {
                const double xk         = r(k,j);
                const double xi         = r(i,j);
                sum += (xk - xi) / R(i,k) * m_jastrowGradient(i, k);
            }
            for (int i = k + 1; i < m_numberOfElectrons; i++) {
                const double xk         = r(k,j);
                const double xi         = r(i,j);
                sum -= (xi - xk) / R(i,k) * m_jastrowGradient(k,i);
            }
            m_quantumForce(k , j) = 2 * slaterGradient(kSpinIndex,j);
            if (m_jastrow) m_quantumForce(k ,j) += 2 * sum;
            m_energyCrossTerm -= 0.5*sum*sum + (slaterGradient(kSpinIndex,j) * sum);
        }
    }
}


void SlaterWithJastrow::evaluateWaveFunctionInitial() {
    if (m_orbitalSet==false) {
        std::cout << "No orbital to use in SlaterWithJastrow wavefunction set." << std::endl;
        exit(1);
    }
    if (m_initializeDone) {
        return;
    }
    m_numberOfSpinUpElectrons   = m_system->getSpinUpElectrons().size();
    m_numberOfSpinDownElectrons = m_system->getSpinDownElectrons().size();
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

    m_spinMatrix             = zeros<mat>(m_numberOfElectrons, m_numberOfElectrons);
    m_electronPositions      = zeros<mat>(m_numberOfElectrons, 3);
    m_interElectronDistances = zeros<mat>(m_numberOfElectrons, m_numberOfElectrons);

    for (int i = 0; i < m_numberOfElectrons; i++) {
        Electron* iElectron = m_system->getElectrons().at(i);
        double xi = iElectron->getPosition().at(0);
        double yi = iElectron->getPosition().at(1);
        double zi = iElectron->getPosition().at(2);
        m_electronPositions(i,0) = xi;
        m_electronPositions(i,1) = yi;
        m_electronPositions(i,2) = zi;

        for (int j = 0; j < m_numberOfElectrons; j++) {
            Electron* jElectron = m_system->getElectrons().at(j);


            if (i==j) {
                m_spinMatrix(i,j) = nan("");
                m_interElectronDistances(i,i) = sqrt(xi*xi + yi*yi + zi*zi);
            } else {
                m_spinMatrix(i,j) = spinCoefficient(iElectron, jElectron);
                m_interElectronDistances(i,j) = computeInterEletronDistance(iElectron, jElectron);
            }
        }
    }
    m_electronPositionsOld      = m_electronPositions;
    m_interElectronDistancesOld = m_interElectronDistances;

    m_slaterUp   = zeros<mat>(eUp,   eUp);
    m_slaterDown = zeros<mat>(eDown, eDown);

    vector<Electron*> spinUpElectrons   = m_system->getSpinUpElectrons();
    vector<Electron*> spinDownElectrons = m_system->getSpinDownElectrons();

    for (int electron = 0; electron < eUp; electron++) {
        const double x = spinUpElectrons.at(electron)->getPosition().at(0);
        const double y = spinUpElectrons.at(electron)->getPosition().at(1);
        const double z = spinUpElectrons.at(electron)->getPosition().at(2);

        for (int basis = 0; basis < eUp; basis++) {
            m_slaterUp(electron, basis) = m_orbital->evaluate(x,y,z,basis,1);
        }
    }
    for (int electron = 0; electron < eDown; electron++) {
        const double x = spinDownElectrons.at(electron)->getPosition().at(0);
        const double y = spinDownElectrons.at(electron)->getPosition().at(1);
        const double z = spinDownElectrons.at(electron)->getPosition().at(2);

        for (int basis = 0; basis < eDown; basis++) {
            m_slaterDown(electron, basis) = m_orbital->evaluate(x,y,z,basis,0);
        }
    }

    m_slaterUp      = m_slaterUp.i();
    m_slaterDown    = m_slaterDown.i();
    m_slaterUpOld   = m_slaterUp;
    m_slaterDownOld = m_slaterDown;

    m_correlationMatrix    = zeros<mat>(m_numberOfElectrons, m_numberOfElectrons);
    m_correlationMatrixOld = zeros<mat>(m_numberOfElectrons, m_numberOfElectrons);
    fillCorrelationMatrix();
    m_correlationMatrixOld = m_correlationMatrix;

    m_slaterGradientUp   = zeros<mat>(eUp,   3);
    m_slaterGradientDown = zeros<mat>(eDown, 3);
    m_jastrowGradient       = zeros<mat>(m_numberOfElectrons,m_numberOfElectrons);
    m_jastrowLaplacianTerms = zeros<mat>(m_numberOfElectrons,m_numberOfElectrons);

    for (int electron = 0; electron < m_numberOfElectrons; electron++) {
        updateSlaterGradient(1.0, electron);
        computeSlaterLaplacian(electron);   // This is called N-2 too many times,
        // but that shouldnt affect runtime
        // in any substantial way.
        if (m_jastrow) {
            updateJastrowGradient(electron);
            updateJastrowLaplacianTerms(electron);
            computeJastrowLaplacian();
        }
    }
    m_jastrowGradientOld        = m_jastrowGradient;
    m_jastrowLaplacianTermsOld  = m_jastrowLaplacianTerms;
    m_slaterGradientUpOld       = m_slaterGradientUp;
    m_slaterGradientDownOld     = m_slaterGradientDown;

    m_quantumForce = zeros<mat>(m_numberOfElectrons,3);
    computeQuantumForce();
    m_quantumForceOld = m_quantumForce;

    m_initializeDone = true;
}

void SlaterWithJastrow::passProposedChangeToWaveFunction(int electronChanged, int dimensionChanged) {
    m_changedElectron   = electronChanged;
    m_changedDimension  = dimensionChanged;
    m_spinChanged       = m_system->getElectrons().at(electronChanged)->getSpin();
}

void SlaterWithJastrow::updateWaveFunctionAfterAcceptedStep() {
    m_electronPositionsOld      = m_electronPositions;
    m_quantumForceOld           = m_quantumForce;
    m_correlationMatrixOld      = m_correlationMatrix;
    m_jastrowGradientOld        = m_jastrowGradient;
    m_jastrowLaplacianTermsOld  = m_jastrowLaplacianTerms;
    mat& slaterGradient         = (m_spinChanged==1 ? m_slaterGradientUp    : m_slaterGradientDown);
    mat& slaterGradientOld      = (m_spinChanged==1 ? m_slaterGradientUpOld : m_slaterGradientDownOld);
    slaterGradientOld           = slaterGradient;

    updateSlaterInverse();
    if (m_spinChanged==1) {
        m_slaterUpOld = m_slaterUp;
    } else {
        m_slaterDownOld = m_slaterDown;
    }
    m_interElectronDistancesOld = m_interElectronDistances;
}

void SlaterWithJastrow::updateWaveFunctionAfterRejectedStep() {
    m_electronPositions         = m_electronPositionsOld;
    m_correlationMatrix         = m_correlationMatrixOld;
    m_slaterGradientUp          = m_slaterGradientUpOld;
    m_slaterGradientDown        = m_slaterGradientDownOld;
    m_jastrowGradient           = m_jastrowGradientOld;
    m_jastrowLaplacianTerms     = m_jastrowLaplacianTermsOld;
    m_quantumForce              = m_quantumForceOld;
    m_interElectronDistances    = m_interElectronDistancesOld;
}

double SlaterWithJastrow::computeWaveFunctionRatio(int changedElectronIndex) {
    m_changedElectron = changedElectronIndex;
    updateElectronDistanceMatrices();
    updateCorrelationsMatrix();
    updateJastrowGradient(m_changedElectron);
    computeJastrowRatio();
    computeSlaterRatio();
    updateSlaterGradient(m_Rsd, m_changedElectron);
    computeQuantumForce();

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





