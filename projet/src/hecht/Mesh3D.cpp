// ====================================
// Frédéric Hecht - Sorbonne Université
// ====================================

#include <fstream>
#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include "hecht/Mesh3D.hpp"

using namespace std;

Mesh3d::Mesh3d(string filename) : filename(filename) {

    std::ifstream f(filename);
    assert(f);

    int i, I5[5];
    string str;

    while (!f.eof()) {

        f >> str;

        if (str == "Vertices") {

            f >> vertexCount;
            assert(!this->vertices);
//            cout << "  -- Nb of Vertex " << vertexCount << endl;
            this->vertices = new Vertex3[vertexCount];

            for (i = 0; i < vertexCount; i++) {

                f >> vertices[i];
                vertices[i].index = i;
                assert(f.good());

                if (i == 0) {
                    minimum = vertices[i];
                    maximum = vertices[i];

                } else {
                    minimum.x = min(minimum.x, vertices[i].x);
                    minimum.y = min(minimum.y, vertices[i].y);
                    minimum.z = min(minimum.z, vertices[i].z);
                    maximum.x = max(maximum.x, vertices[i].x);
                    maximum.y = max(maximum.y, vertices[i].y);
                    maximum.z = max(maximum.z, vertices[i].z);
                }
            }

        } else if (str == "Tetrahedra") {

            f >> elementCount;
            assert(this->vertices && !this->elements);
            this->elements = new Simplex3[elementCount];
            mes = 0;
            assert(this->elements);

//            cout << "  -- Nb of Elements " << elementCount << endl;

            for (int i = 0; i < elementCount; i++) {

                for (int k = 0; k < 5; ++k) {
                    f >> I5[k];
                }

                this->elements[i].build(this->vertices, I5, -1);
                mes += this->elements[i].volume;
            }

        } else if (str == "Triangles") {

            int triangleCount;
            f >> triangleCount;

            for (i = 0; i < triangleCount; i++) {

                int vertexIndices[3];
                int label;

                f >> vertexIndices[0] >> vertexIndices[1] >> vertexIndices[2] >> label;

                for (int i = 0; i < 3; i++) {
                    vertices[vertexIndices[i] - 1].boundaryLabels.push_back(label);
                }
            }

        } else if (str[0] == '#') {

            int c;

            while ((c = f.get()) != '\n' && c != EOF) {
                ;
            }
        }
    }

    assert (elementCount >= 0 && vertexCount > 0);
//    std::cout << " End read " << vertexCount << " " << elementCount << " mes =" << mes << std::endl;
}
