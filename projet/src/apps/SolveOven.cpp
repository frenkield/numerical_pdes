#include "solver/HelmholtzSolver3D.hpp"
#include "solver/HeatSolver3D.hpp"
#include "OvenParameters.cpp"

using namespace std;

int main(int argc, const char **argv) {

    // on commence avec l'equation de helmholtz...

    string meshFilename = objectInOven ? "four_3d.mesh" : "four_vide_3d.mesh";
    HelmholtzSolver3D helmholtzSolver("maillages/" + meshFilename);

    helmholtzSolver.reflectionConditions = helmholtzReflectionConditions;

    // il faut calibrer les ondes sur l'axe x
    int ovenLength = helmholtzSolver.mesh->maximum.x - helmholtzSolver.mesh->minimum.x;
    helmholtzSolver.angularFrequency = wavelengthCount * 2 * M_PI / ovenLength;

    helmholtzSolver.computePermittivity = computePermittivity;
    helmholtzSolver.computeBoundaryValue = computeBoundaryValuePyramid;

    cout << endl << "Calcul de la solution de l'equation de helmholtz" << endl;
    helmholtzSolver.solve();
    helmholtzSolver.displayStatistics();
    helmholtzSolver.writeDataForPlots();

    vector<complex<double>>& helmholtzResult = helmholtzSolver.u;

    // ===================================================================

    // et puis l'equation de poisson...

    HeatSolver3D heatSolver(helmholtzSolver.mesh);
    heatSolver.computeThermalConductivity = computeThermalConductivity;

    vector<double> helmholtzMagnitude(helmholtzSolver.mesh->vertexCount, 0);

    // calcul de u * u_bar
    transform(helmholtzResult.begin(), helmholtzResult.end(), helmholtzMagnitude.begin(),
              [](complex<double> d) -> double {return norm(d) / 500;});

    heatSolver.buildRightHandSide(helmholtzMagnitude);

    cout << endl << "Calcul de la solution de l'equation de poisson" << endl;
    heatSolver.solve();
    heatSolver.displayStatistics();
    heatSolver.writeDataForPlots();

    // on affiche les résultats avec FreeFEM...

    if (objectInOven) {
        heatSolver.runFreefemScript("view_oven_solution.edp", "-objectInOven 1");
    } else {
        heatSolver.runFreefemScript("view_oven_solution.edp", "-objectInOven 0");
    }
}
