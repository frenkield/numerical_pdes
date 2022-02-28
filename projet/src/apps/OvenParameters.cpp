#include "solver/HelmholtzSolver3D.hpp"

using namespace std;

// ======================= toutes les constantes =========================

const bool objectInOven = true;
const bool helmholtzReflectionConditions = false;

const double wavelengthCount = 1.85;

const complex<double> permittivityAir = {-1.0, 0.05};
const complex<double> inversePermittivityAir = 1.0 / permittivityAir;
const complex<double> inversePermittivityObject = inversePermittivityAir * 0.25;

const double thermalConductivityAir = 0.0262;
const double thermalConductivityObject = 0.6;

// =================== calculateurs des non-constantes =====================

bool inObject(Simplex3 simplex) {
    return objectInOven && simplex.label == 1;
}

complex<double> computePermittivity(Simplex3& simplex) {

    if (inObject(simplex)) {
        return inversePermittivityObject;
    }

    return inversePermittivityAir;
}

double computeThermalConductivity(Simplex3& simplex) {

    if (inObject(simplex)) {
        return thermalConductivityObject;
    }

    return thermalConductivityAir;
}

/**
 * Le calcul d'une valeur sur le bord s'effectue avec des points
 * normalisés. C'est-à-dire, l'argument normalizedVertex contient
 * les coordonnées normalisées de sorte que x, y, et z aient des
 * valuers entre 0 et 1.
 */
double computeBoundaryValuePyramid(R3& normalizedVertex) {

    double slope = 4;

    double boundaryValue = min({
        slope * normalizedVertex.z,
        slope * normalizedVertex.y,
        slope * (1 - normalizedVertex.z),
        slope * (1 - normalizedVertex.y),
        1.0
    });

    return boundaryValue;
}
