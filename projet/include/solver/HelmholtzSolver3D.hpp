#ifndef SOLVER_HELMHOLTZSOLVER3D_H
#define SOLVER_HELMHOLTZSOLVER3D_H

#include "Solver3D.hpp"

typedef complex<double> TypeScalar;

class HelmholtzSolver3D : public Solver3D<TypeScalar> {

public:

    HelmholtzSolver3D(string meshFilename);
    double angularFrequency = 1.0;
    bool reflectionConditions = false;

    void configureOperations();
    void applyPenalization();

    function<complex<double>(Simplex3&)> computePermittivity = [](Simplex3 simplex) {
        return (complex<double>){1, 1};
    };

    function<double(R3&)> computeBoundaryValue = [](R3& vertex) {
        return 0.0;
    };
};

#endif