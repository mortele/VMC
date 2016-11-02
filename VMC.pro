TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -L/usr/local/lib -larmadillo -llapack -lblas
INCLUDEPATH += /usr/local/include

SOURCES += main.cpp \
    Hamiltonians/hamiltonian.cpp \
    system.cpp \
    electron.cpp \
    atom.cpp

HEADERS += \
    Hamiltonians/hamiltonian.h \
    system.h \
    electron.h \
    atom.h
