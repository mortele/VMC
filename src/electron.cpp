#include "electron.h"

Electron::Electron(arma::vec position, int spin) {
    m_spin      = spin;
    m_position  = position;
}

void Electron::setPosition(arma::vec position) {
    m_position = position;
}

void Electron::setPosition(double position, int dimension) {
    m_position(dimension) = position;
}

void Electron::adjustPosition(double adjustment, int dimension) {
    m_position(dimension) += adjustment;
}
