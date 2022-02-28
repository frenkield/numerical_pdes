// ====================================
// Frédéric Hecht - Sorbonne Université
// ====================================

#ifndef HECHT_MESH3D_HPP
#define HECHT_MESH3D_HPP

#include <cassert>
#include "hecht/R3.hpp"
#include "hecht/Label.hpp"
#include "Simplex3.hpp"

class Mesh3d {

public:

    int vertexCount = 0, elementCount = 0;
    Vertex3 *vertices = 0;
    Simplex3 *elements = 0;
    double mes = 0;
    R3 minimum;
    R3 maximum;
    string filename;

    Mesh3d(std::string filename);

    ~Mesh3d() {
        delete[] vertices;
        delete[] elements;
    }

    // destuctor => careful with copie operator
    // no copy operator
    // chech index number
    int CheckV(int i) const {
        assert(i >= 0 && i < vertexCount);
        return i;
    }

    int CheckT(int i) const {
        assert(i >= 0 && i < elementCount);
        return i;
    }

    int operator()(const Vertex3 &vv) const { return CheckV(&vv - vertices); }
    int operator()(const Simplex3 &tt) const { return CheckT(&tt - elements); }

    int operator()(const Vertex3 *vv) const { return CheckV(vv - vertices); }  // (1)
    int operator()(const Simplex3 *tt) const { return CheckT(tt - elements); }

    Simplex3 &operator[](int k) { return elements[CheckT(k)]; }
    const Simplex3 &operator[](int k) const { return elements[CheckT(k)]; }

    int operator()(int k, int i) const { return operator()(elements[k].v[i]); }// call (1)

private:
    Mesh3d(const Mesh3d &);
    Mesh3d &operator=(const Mesh3d &);
};

#endif
