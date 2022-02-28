#ifndef SOLVER_SOLVER3D_H
#define SOLVER_SOLVER3D_H

#include <vector>
#include <complex>
#include <chrono>
#include "hecht/EF3d.hpp"
#include "hecht/Mesh3D.hpp"
#include "Timer.hpp"

using namespace std;
using namespace std::chrono;

template<class TypeScalar>
struct Operation {
    function<TypeScalar(Simplex3&)> computeWeight;
    int operatorU, operatorV;
};

/**
 * Résoudre par des éléments finis une equation differentielle partielle.
 */
template<class TypeScalar>
class Solver3D {

    typedef HashMatrix<int, TypeScalar> HashMatrixType;
    typedef Operation<TypeScalar> OperationType;

public:

    string filenamePrefix = "solver3d";
    shared_ptr<Mesh3d> mesh;
    unique_ptr<HashMatrixType> A;
    vector<OperationType> operations;
    vector<TypeScalar> u;
    vector<TypeScalar> b;

    Solver3D(string filenamePrefix, string meshFilename) : filenamePrefix(filenamePrefix) {

        cout << "Lecture de " << meshFilename << endl;
        mesh = unique_ptr<Mesh3d>(new Mesh3d(meshFilename));

        A = unique_ptr<HashMatrixType>(new HashMatrixType(mesh->vertexCount));

        u.resize(mesh->vertexCount, 0);
        b.resize(mesh->vertexCount, 0);
    }

    Solver3D(string filenamePrefix, shared_ptr<Mesh3d> existingMesh) : filenamePrefix(filenamePrefix) {
        this->mesh = existingMesh;
        A = unique_ptr<HashMatrixType>(new HashMatrixType(mesh->vertexCount));
        u.resize(mesh->vertexCount, 0);
        b.resize(mesh->vertexCount, 0);
    }

    virtual void applyPenalization() {}
    virtual void configureOperations() {}

    void solve() {

        Timer timer;

        configureOperations();

        timer.reset();
        buildMatrix();
        applyPenalization();
        timer.logElapsedTime("Temps assemblage matrice");

        double eps = 1e-8;

        timer.reset();
        SparseLinearSolver<int, TypeScalar> solver(*A, "UMFPACK", "eps", eps, nullptr);
        solver.solve(&u[0], &b[0]);
        timer.logElapsedTime("Temps solve umfpack");
    }

    void addOperation(OperationType operation) {
        operations.push_back(operation);
    }

    static TypeScalar weightOne(Simplex3 simplex) {
        return 1;
    }

    void addOperation(int uOperation, int vOperation) {
        addOperation(weightOne, uOperation, vOperation);
    }

    void addOperation(function<TypeScalar(Simplex3&)> weight, int uOperation, int vOperation) {
        addOperation({weight, uOperation, vOperation});
    }

    void buildMatrix(Simplex3& simplex, const vector<double>& D, const vector<double>& W) {

        int npq = D.size();
        int vertexCount = simplex.vertexCount;

        for (int ip = 0; ip < vertexCount; ++ip) {

            for (int jp = 0; jp < vertexCount; ++jp) {

                Vertex3 vertexI = simplex[ip];
                Vertex3 vertexJ = simplex[jp];

                TypeScalar aKij = 0;

                for (int op = 0; op < operations.size(); ++op) {

                    OperationType operation = operations[op];
                    TypeScalar weight = operation.computeWeight(simplex);

                    for (int p = 0; p < npq; ++p) {

                        int k = p * vertexCount * vertexCount;

                        int indexV = k + ip * vertexCount + operation.operatorV;
                        int indexU = k + jp * vertexCount + operation.operatorU;

                        aKij += D[p] * weight * W[indexV] * W[indexU];
                    }
                }

                (*A)(vertexI.index, vertexJ.index) += aKij;
            }
        }
    }

    void buildMatrix() {

        const QF3& qf = QuadratureFormular_Tet_2;
        vector<double> Wh, Dh;
        SetWP1(qf, Wh, Dh);

        vector<double> WK(Wh.size()), DK(Dh.size());

        for (int k = 0; k < mesh->elementCount; ++k) {
            Simplex3& simplex = (*mesh)[k];
            SetWK(simplex, Wh, WK);
            SetDK(simplex, Dh, DK);
            buildMatrix(simplex, DK, WK);
        }
    }

    void runFreefemScript(string scriptFilename, string arguments = "") {
        string command = "FreeFem++ freefem/" + scriptFilename + " -ns " + arguments;
        cout << "Execution de freefem : " << command << endl;
        system(command.c_str());
    }

    void writeDataForPlots() {

        cout << "Écriture des données pour ParaView et FreeFem" << endl;
        ofstream fichierFreeFem("visualisation/freefem_" + filenamePrefix + ".txt");

        for (int k = 0; k < mesh->elementCount; k++) {

            Simplex3& simplex = mesh->elements[k];

            for (int ip = 0; ip < 4; ip++) {
                complex<double> value = u[simplex.v[ip]->index];
                fichierFreeFem << value << " ";
            }

            fichierFreeFem << endl;
        }

        // ================================================================

        ofstream fichierParaview("visualisation/" + filenamePrefix + ".csv");
        fichierParaview << "x,y,z,real,imag,mag" << endl;

        for (int k = 0; k < mesh->vertexCount; k++) {

            Vertex3& vertex = mesh->vertices[k];

            complex<double> value = u[vertex.index];
            double magnitude = norm(value);

            fichierParaview << vertex.x << "," << vertex.y << "," << vertex.z << ",";
            fichierParaview << value.real() << "," << value.imag() << "," << magnitude << endl;
        }
    }

    void displayStatistics() {

        double maxMagnitude = 0;
        double averageMagnitude = 0;

        for (int k = 0; k < mesh->vertexCount; k++) {

            Vertex3& vertex = mesh->vertices[k];

            complex<double> value = u[vertex.index];
            double magnitude = norm(value);

            averageMagnitude += magnitude;
            maxMagnitude = max(maxMagnitude, magnitude);

            if (vertex.x >= -0.1 && vertex.x <= 0.1
                && vertex.y >= -0.1 && vertex.y <= 0.1
                && vertex.z >= -0.1 && vertex.z <= 0.1) {
                cout << "Magnitude au centre = " << magnitude << endl;
            }
        }

        averageMagnitude /= mesh->vertexCount;
        cout << "Max magnitude = " << maxMagnitude << ", ";
        cout << "Magnitude moyenne = " << averageMagnitude << endl;
    }

    static string getMeshFilename(int argc, const char *const *argv) {

        string meshFilename = "maillages/four_3d.mesh";

        if (argc > 1) {
            meshFilename = argv[1];
        }

        return meshFilename;
    }
};

#endif