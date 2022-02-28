#include "solver/HelmholtzSolver3D.hpp"

using namespace std;

HelmholtzSolver3D::HelmholtzSolver3D(string meshFilename)
    : Solver3D("helmholtz", meshFilename) {}

void HelmholtzSolver3D::applyPenalization() {

    double penalizationFactor = 1e30;
    R3 size = mesh->maximum - mesh->minimum;

    for (int k = 0; k < mesh->vertexCount; k++) {

        Vertex3& vertex = mesh->vertices[k];

        if (!reflectionConditions) {
            // dirichlet homogene sur (presque) toutes les parois du four
            if (vertex.isOnBoundary() && !vertex.isOnBoundary(7)) {
                (*A)(k, k) = penalizationFactor;
            }
        }

        // applique les conditions aux limites g(x) sur la paroi x = -20
        if (vertex.isOnBoundary(1)) {

            (*A)(k, k) = penalizationFactor;
            R3 normalizedVertex = vertex - mesh->minimum;

            normalizedVertex.x /= size.x;
            normalizedVertex.y /= size.y;
            normalizedVertex.z /= size.z;

            double boundaryValue = computeBoundaryValue(normalizedVertex);

            b[k] = penalizationFactor * boundaryValue;
            u[k] = boundaryValue;
        }
    }
}

/**
 * Configurer l'operations pour la forme variationnelle de l'EDP.
 *
 * L(u,v) = (omega^2*u*v) + (1/epsilon * (du/dx*dv/dx + du/dy*dv/dy + du/dz*dv/dz))
 */
void HelmholtzSolver3D::configureOperations() {

    complex<double> angularFrequencySquared = angularFrequency * angularFrequency;

    auto getAngularFrequencySquared = [angularFrequencySquared](Simplex3 simplex) {
        return angularFrequencySquared;
    };

    // w^2 * u * v (on suppose mu = 1)
    addOperation(getAngularFrequencySquared, 0, 0);

    // 1/epsilon * du/dx * dv/dx
    addOperation(computePermittivity, 1, 1);

    // 1/epsilon * du/dy * dv/dy
    addOperation(computePermittivity, 2, 2);

    // 1/epsilon * du/dz * dv/dz
    addOperation(computePermittivity, 3, 3);
}
