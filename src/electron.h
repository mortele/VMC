#pragma once
#include <armadillo>
#include <iostream>

class Electron {
private:
    int       m_spin;
    int       m_index;
    int       m_spinIndex;
    arma::vec m_position;

public:
    Electron(arma::vec position, int spin=1);

    void setPosition    (arma::vec position);
    void setPosition    (double position, int dimension);
    void adjustPosition (double adjustment, int dimension);

    int         getSpin()       { return m_spin;        }
    arma::vec   getPosition()   { return m_position;    }

    void  setIndex    (int i)   { m_index       = i;    }
    void  setSpinIndex(int i)   { m_spinIndex   = i;    }
    int   getSpinIndex()        { return m_spinIndex;   }
    int   getIndex()            { return m_index;       }
};

