TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -L/usr/local/lib -larmadillo -llapack -lblas
INCLUDEPATH += /usr/local/include

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
    WaveFunctions/directevaluationslaterwithjastrow.cpp

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
    WaveFunctions/directevaluationslaterwithjastrow.h
