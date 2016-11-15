TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -L/usr/local/lib -larmadillo -llapack -lblas
INCLUDEPATH += /usr/local/include

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

SOURCES += main.cpp \
    hamiltonian.cpp \
    system.cpp \
    electron.cpp \
    sampler.cpp \
    WaveFunctions/wavefunction.cpp \
    WaveFunctions/hydrogenwavefunction.cpp \
    Cores/atom.cpp \
    Cores/core.cpp \
    metropolis.cpp \
    WaveFunctions/heliumwavefunction.cpp \
    unittest.cpp \
    WaveFunctions/heliumwithjastrow.cpp \
    WaveFunctions/directevaluationslater.cpp \
    WaveFunctions/directevaluationslaterwithjastrow.cpp \
    hartreefockbasisparser.cpp \
    WaveFunctions/Orbitals/primitivegaussian.cpp \
    WaveFunctions/Orbitals/contractedgaussian.cpp \
    WaveFunctions/gaussianslater.cpp

HEADERS += \
    hamiltonian.h \
    system.h \
    electron.h \
    sampler.h \
    WaveFunctions/wavefunction.h \
    WaveFunctions/hydrogenwavefunction.h \
    Cores/atom.h \
    Cores/core.h \
    metropolis.h \
    RandomNumberGenerator/random.h \
    WaveFunctions/heliumwavefunction.h \
    unittest.h \
    WaveFunctions/heliumwithjastrow.h \
    WaveFunctions/directevaluationslater.h \
    WaveFunctions/directevaluationslaterwithjastrow.h \
    hartreefockbasisparser.h \
    WaveFunctions/Orbitals/primitivegaussian.h \
    WaveFunctions/Orbitals/contractedgaussian.h \
    WaveFunctions/gaussianslater.h
