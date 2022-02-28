// ====================================
// Frédéric Hecht - Sorbonne Université
// ====================================

#ifndef HECHT_EF3D_VERTEX_HPP
#define HECHT_EF3D_VERTEX_HPP

#include <vector>
#include "hecht/R3.hpp"
#include "hecht/Label.hpp"

class Vertex3 : public R3, public Label {

public:

    Vertex3() : R3() {}
    Vertex3(double x, double y, double z) : R3(x, y, z) {}
    int index = -1;
    vector<int> boundaryLabels;

    bool isOnBoundary(int boundaryLabel) {
        return find(boundaryLabels.begin(), boundaryLabels.end(), boundaryLabel) != boundaryLabels.end();
    }

    bool isOnBoundary() {
        return boundaryLabels.size() > 0;
    }
};

inline std::ostream &operator<<(std::ostream &f, const Vertex3 &P) {
    return f << P.x << ' ' << P.y << ' ' << P.z << ' ' << P.lab << ' ';
}

inline std::istream &operator>>(std::istream &f, Vertex3 &P) { return f >> P.x >> P.y >> P.z >> P.lab; }

#endif
