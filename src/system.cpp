#include "system.h"
#include "electron.h"
#include "hamiltonian.h"
#include "metropolis.h"
#include "sampler.h"
#include "optimizer.h"
#include "WaveFunctions/wavefunction.h"
#include "Cores/core.h"

System::System() {
    m_numberOfElectrons         = 0;
    m_numberOfSpinUpElectrons   = 0;
    m_numberOfSpinDownElectrons = 0;
    m_numberOfDimensions        = 3;
    m_sampler                   = new Sampler(this);
    m_metropolis                = new Metropolis(this);
    m_hamiltonian               = new Hamiltonian(this);
    m_optimizer                 = new Optimizer(this);
}

void System::setup() {
    m_hamiltonian->setup();
    m_hamiltonian->setElectronInteraction(m_interactingElectrons);
    m_metropolis->setup();
    m_sampler->setup();
    setupSpinElectronArrays();
    m_waveFunction->setup();
    m_optimizer->setup();
}

void System::setupSpinElectronArrays() {
    for (int electron = 0; electron < m_numberOfElectrons; electron++) {
        int spin = m_electrons.at(electron)->getSpin();

        // Spin up electron.
        if (spin == 1) {
            m_spinUpElectrons.push_back(m_electrons.at(electron));

        // Spin down electron.
        } else if (spin == 0) {
            m_spinDownElectrons.push_back(m_electrons.at(electron));

        } else {
            std::cout << "Something is terribly wrong in System::setupSpinElectronArrays()." << std::endl;
        }
    }
    m_numberOfSpinUpElectrons   = m_spinUpElectrons.size();
    m_numberOfSpinDownElectrons = m_spinDownElectrons.size();
}

double System::runMetropolis(int steps) {
    setup();
    double E = m_metropolis->runSteps(steps);
    return E;
}

double System::runMetropolisSilent(int steps) {
    setup();
    double E = m_metropolis->runStepsSilent(steps);
    return E;
}

double System::optimizeBeta(double beta,
                            double tollerance,
                            int    maxIterations,
                            int    cycles) {
    m_optimizer->optimizeBeta(beta,tollerance,maxIterations,cycles);
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setElectronInteraction(bool interacting) {
    m_interactingElectrons = interacting;
}

void System::setStepLength(double stepLength) {
    m_metropolis->setStepLength(stepLength);
}

void System::setImportanceSampling(bool importanceSampling) {
    m_metropolis->setImportanceSampling(importanceSampling);
}

void System::setOrbital(Orbital* orbital) {
    m_waveFunction->setOrbital(orbital);
}

void System::adjustPositionOfElectron(int       electron,
                                      int       dimension,
                                      double    adjustment) {
    m_electrons.at(electron)->adjustPosition(adjustment, dimension);
}

void System::addCore(Core* core) {
    m_cores.push_back(core);
    for (Electron* electron : core->getElectrons()) {
        m_electrons.push_back(electron);
        m_numberOfElectrons += 1;
    }
}
