#include "solver/HeatSolver3D.hpp"
#include "solver/TrivialSolver3D.hpp"

using namespace std;

HeatSolver3D::HeatSolver3D(string meshFilename) : Solver3D("heat", meshFilename) {
}

HeatSolver3D::HeatSolver3D(shared_ptr<Mesh3d> mesh) : Solver3D("heat", mesh) {
}

void HeatSolver3D::buildRightHandSide(vector<double>& rightHandSide) {

    TrivialSolver3D rightHandSideSolver(mesh);
    rightHandSideSolver.buildMatrix();

    rightHandSideSolver.A->addMatMul(&rightHandSide[0], &b[0]);
}

void HeatSolver3D::applyPenalization() {

    double tgv = 1e30;
    R3 size = mesh->maximum - mesh->minimum;

    for (int k = 0; k < mesh->vertexCount; k++) {

        Vertex3& vertex = mesh->vertices[k];

        // il faut ignorer la viande ici
        if (vertex.isOnBoundary() && !vertex.isOnBoundary(7)) {

            (*A)(k, k) = tgv;

            double boundaryValue = 0;
            b[k] = tgv * boundaryValue;
            u[k] = boundaryValue;
        }
    }
}

void HeatSolver3D::configureOperations() {

    // u*v + weight * (du/dx*dv/dx + du/dy*dv/dy + du/dz*dv/dz)

    // weight * du/dx * dv/dx
    addOperation(computeThermalConductivity, 1, 1);

    // weight * du/dy * dv/dy
    addOperation(computeThermalConductivity, 2, 2);

    // weight * du/dz * dv/dz
    addOperation(computeThermalConductivity, 3, 3);
}
