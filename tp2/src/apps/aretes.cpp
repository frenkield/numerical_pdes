#include <iostream>
#include <vector>
#include <cmath>
#include "../EF2d-base.hpp"

using namespace std;


int main(int argc, const char **argv) {

    Mesh2d mesh1("mesh/mesh1.msh");

    //~ for (int i = 0; i < mesh1.nv; i++) {
        //~ cout << mesh1.v[i] << endl;
    //~ }

    for(int k = 0; k < mesh1.nt; ++k) {

        cout << mesh1(k, 0) << " " << mesh1(k, 1) << " " <<
                mesh1(k, 2) << endl;

        //~ for(int ip = 0; ip <= 3; ++ip) {
            //~ int i3 = ip % 3;
            //~ int i = mesh1(k, i3);
            //~ cout << i << endl;
        //~ }
    }


    return 0;
}
