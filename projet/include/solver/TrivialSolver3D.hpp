#ifndef SOLVER_TRIVIALSOLVER3D_H
#define SOLVER_TRIVIALSOLVER3D_H

#include "solver/Solver3D.hpp"

typedef double TypeScalar;

/**
 * Ce n'est pas vraiment un solveur mais ça facilite des choses.
 *
 * On utilise cette classe pour générer la matrice de masse afin de calculer le
 * membre de droite.
 */
class TrivialSolver3D : public Solver3D<TypeScalar> {

public:

    TrivialSolver3D(string meshFilename) : Solver3D("trivial", meshFilename) {
        addOperation(0, 0);
    }

    TrivialSolver3D(shared_ptr<Mesh3d> mesh) : Solver3D("trivial", mesh) {
        addOperation(0, 0);
    }
};

#endif