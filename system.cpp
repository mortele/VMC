#include "system.h"
#include "electron.h"
#include "hamiltonian.h"
#include "metropolis.h"
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "Cores/core.h"
#include "RandomNumberGenerator/random.h"

System::System() {
    m_numberOfElectrons         = 0;
    m_numberOfSpinUpElectrons   = 0;
    m_numberOfSpinDownElectrons = 0;
    m_numberOfDimensions        = 3;
    m_sampler                   = new Sampler(this);
    m_metropolis                = new Metropolis(this);
    m_hamiltonian               = new Hamiltonian(this);
    Random::seed(42413213);
}

void System::setup() {
    m_waveFunction->setup();
    m_hamiltonian->setup();
    m_hamiltonian->setElectronInteraction(m_interactingElectrons);
    m_metropolis->setup();
    m_sampler->setup();
}

void System::runMetropolis(int steps) {
    setup();
    m_metropolis->runSteps(steps);
}

void System::runMetropolisSilent(int steps) {
    setup();
    m_metropolis->runStepsSilent(steps);
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
