#include "core.h"
#include "system.h"


Core::Core(System* system, arma::vec position) {
    m_position  = position;
    m_system    = system;
}

