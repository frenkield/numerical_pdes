// ====================================
// Frédéric Hecht - Sorbonne Université
// ====================================

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <vector>
#include "GC.hpp"

using namespace std;


class MatFullRowMajor : public MatVirt {
public:
    int n, m; //
    double *a;

    MatFullRowMajor(int nn, int mm, double *v = 0) : MatVirt(nn, mm), n(nn), m(mm), a(new double[n * m]) {
        if (v) for (int k = 0; k < n * m; ++k) a[k] = v[k];
        else std::fill(a, a + n * m, 0.);
    }

    //  Pb de destruction car "v" n'est pas forcement alloue avec new
//    MatFullRowMajor(int nn,int mm,double *v) : n(nn),m(mm),a(v) {}
    double &operator()(int i, int j) { return a[m * i + j]; }

    const double &operator()(int i, int j) const { return a[m * i + j]; }

    ~MatFullRowMajor() { delete[] a; }

    MatFullRowMajor(const MatFullRowMajor &A) : MatVirt(A.n, A.m), n(A.n), m(A.m), a(new double[n * m]) {
        operator=(A);
    }

    MatFullRowMajor &operator=(const MatFullRowMajor &A) {
        assert(n == A.n && m == A.m);
        int nm = n * m;
        for (int k = 0; k < nm; ++k)
            a[k] = A.a[k];
        return *this;
    }

    double *addmatmul(double *x, double *Ax) const {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                Ax[i] += a[i * m + j] * x[j];
        return Ax;
    }
};


//int main(int argc, const char ** argv)
//{
//    int n = 4;
//    double  acoef[] = {10.,1.,2.,4.,
//        1.,11.,0.,3.,
//        2.,0.,12.,-1.,
//        4.,3.,-1.,13.
//    } ;
//    MatFullRowMajor A(n,n,acoef);
//    vector<double> x0(n),b(n),x(n,0.);
//    for(int i=0; i< n; ++i) x0[i]=i+1;
//    for( int i=0; i< n; ++i)
//    {
//        for( int j=0; j< n; ++j)
//            cout << A(i,j) << " ";
//            cout << endl;
//    }
//    MatIdentite Id(n);
//    ProduitMatVec(A,&x0[0],&b[0]);
//    GradientConjugue(A,Id,&b[0],&x[0],n,1e-8,10);
//    cout << "b=" << endl;
//    for( int i=0; i<n; ++i) cout <<b[i] << "\t" << x[i] << "\t" <<  endl;
//   return 0;
//}
