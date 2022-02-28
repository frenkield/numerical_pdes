// ====================================
// Frédéric Hecht - Sorbonne Université
// ====================================

#ifndef HECHT_EF3D_SIMPLEX3_HPP
#define HECHT_EF3D_SIMPLEX3_HPP

#include "Vertex3.hpp"
#include <cassert>

class Simplex3 {

public:

    static const int vertexCount = 4;
    Vertex3 *v[vertexCount];
    double volume;
    int label = -1;

    Simplex3() {
        (v[0] = (v[1] = (v[2] = 0)));
    }

    void build(Vertex3 *v0, int *I, int offset = 0) {// I array of vertex number

        for (int i = 0; i < vertexCount; ++i) {
            v[i] = v0 + I[i] + offset;
        }

        label = I[4];

        volume = det(*v[0], *v[1], *v[2], *v[3]) / 6.;
        assert(volume > 0);
    }

    void GradLambdaK(R3 *G) const {
        double K6 = volume * 6;
        G[0] = -(R3(*v[1], *v[2]) ^ R3(*v[1], *v[3])) / K6;
        G[1] = (R3(*v[0], *v[2]) ^ R3(*v[0], *v[3])) / K6;
        G[2] = -(R3(*v[0], *v[1]) ^ R3(*v[0], *v[3])) / K6;
        G[3] = -(G[0] + G[1] + G[2]);
        // G[3] = (R3(*v[0], *v[1])^R3(*v[0], *v[2]))/K6;
    }

    Vertex3 &operator[](int i) {
        assert(i >= 0 && i < vertexCount);
        return *(v[i]);
    }

    const Vertex3 &operator[](int i) const {
        assert(i >= 0 && i < vertexCount);
        return *(v[i]);
    }

    Vertex3 getCenter() const {
        Vertex3 center;
        center.x = (v[0]->x + v[1]->x + v[2]->x + v[3]->x) / 4;
        center.y = (v[0]->y + v[1]->y + v[2]->y + v[3]->y) / 4;
        center.z = (v[0]->z + v[1]->z + v[2]->z + v[3]->z) / 4;
        return center;
    }
};

#endif
