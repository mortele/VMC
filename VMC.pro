TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -L/usr/local/lib -larmadillo -llapack -lblas
INCLUDEPATH += /usr/local/include

QMAKE_CXXFLAGS_RELEASE -= -O2 -fPIC -fpic
QMAKE_CXXFLAGS_RELEASE += -O3 #-ffast-math
QMAKE_CFLAGS_APP -= -fPIC
QMAKE_CXXFLAGS_RELEASE += -fno-pic -O3 # -unroll-aggressive -qopt-prefetch -fast -xCORE-AVX2 -ipo# -prof_use

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
    WaveFunctions/gaussianslater.cpp \
    WaveFunctions/harmonicoscillatorwavefunction.cpp \
    Cores/harmonicoscillator.cpp \
    WaveFunctions/Orbitals/hydrogenorbital.cpp \
    WaveFunctions/slaterwithjastrow.cpp \
    WaveFunctions/Orbitals/orbital.cpp \
    WaveFunctions/Orbitals/gaussianorbital.cpp \
    Optimization/interpolation.cpp \
    Optimization/ap.cpp \
    Optimization/dataanalysis.cpp \
    Optimization/fasttransforms.cpp \
    Optimization/integration.cpp \
    Optimization/linalg.cpp \
    Optimization/optimization.cpp \
    Optimization/solvers.cpp \
    Optimization/specialfunctions.cpp \
    Optimization/statistics.cpp \
    Optimization/alglibinternal.cpp \
    Optimization/alglibmisc.cpp \
    Optimization/diffequations.cpp \
    stofitter.cpp

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
    WaveFunctions/gaussianslater.h \
    Math/exponentialapproximations.h \
    WaveFunctions/harmonicoscillatorwavefunction.h \
    Cores/harmonicoscillator.h \
    WaveFunctions/Orbitals/hydrogenorbital.h \
    WaveFunctions/slaterwithjastrow.h \
    /usr/local/include/armadillo_bits/config.hpp \
    WaveFunctions/Orbitals/orbital.h \
    WaveFunctions/Orbitals/gaussianorbital.h \
    Optimization/stdafx.h \
    Optimization/interpolation.h \
    Optimization/ap.h \
    Optimization/dataanalysis.h \
    Optimization/fasttransforms.h \
    Optimization/integration.h \
    Optimization/linalg.h \
    Optimization/optimization.h \
    Optimization/solvers.h \
    Optimization/specialfunctions.h \
    Optimization/statistics.h \
    Optimization/alglibinternal.h \
    Optimization/alglibmisc.h \
    Optimization/diffequations.h \
    stofitter.h
