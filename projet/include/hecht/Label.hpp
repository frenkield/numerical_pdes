// ====================================
// Frédéric Hecht - Sorbonne Université
// ====================================

#ifndef PROJET_LABEL_HPP
#define PROJET_LABEL_HPP

class Label {

public:

    Label(int l = 0) : lab(l) {}

    int lab;
    int OnGamma() const { return lab; }
};

#endif //PROJET_LABEL_HPP
