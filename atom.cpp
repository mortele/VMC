#include "atom.h"
#include "electron.h"

Atom::Atom(arma::vec position, double charge) {
    m_position  = position;
    m_charge    = charge;
    m_electrons.clear();
    m_electrons.reserve(charge);
}
