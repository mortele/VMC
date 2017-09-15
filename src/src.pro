include(../defaults.pri)

CONFIG   += console c++11
CONFIG   -= app_bundle
CONFIG   -= qt

TEMPLATE = lib

TARGET = VMC


#QMAKE_CXXFLAGS_RELEASE -= -O2 -fPIC -fpic
#QMAKE_CXXFLAGS_RELEASE += -O3 #-ffast-math
#QMAKE_CXXFLAGS_APP -= -fPIC
#QMAKE_CXXFLAGS_RELEASE += -O3 #-fno-pic # -unroll-aggressive -qopt-prefetch -fast -xCORE-AVX2 -ipo# -prof_use

SOURCES += \
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
    stofitter.cpp \
    WaveFunctions/Orbitals/slatertypeorbital.cpp

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
    stofitter.h \
    WaveFunctions/Orbitals/slatertypeorbital.h
