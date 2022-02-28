#ifndef SOLVER_HEATSOLVER3D_H
#define SOLVER_HEATSOLVER3D_H

#include "solver/Solver3D.hpp"

// typedef double TypeScalar;

class HeatSolver3D : public Solver3D<double> {

public:

    HeatSolver3D(string meshFilename);
    HeatSolver3D(shared_ptr<Mesh3d> mesh);

    void configureOperations();
    void applyPenalization();

    void buildRightHandSide(vector<double>& rightHandSide);

    function<double(Simplex3&)> computeThermalConductivity = [](Simplex3& simplex) {
        return 1;
    };

    static double computeWeight(Simplex3 simplex);
};

#endif