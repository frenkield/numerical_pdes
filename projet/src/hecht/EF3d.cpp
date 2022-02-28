// ====================================
// Frédéric Hecht - Sorbonne Université
// ====================================

#include "hecht/EF3d.hpp"

using namespace std;

#include <set>

int debug = 0;

double *pdebugF = 0;
int verbosity = 1;

int SetWP0(const R3 &P, double *w) {
    w[0] = 1;
    w[1] = 0;
    w[2] = 0;
    return 3;
}

int SetWP1(const R3 &P, double *WP1) {
    int k = 0;
    WP1[k++] = 1 - P.x - P.y;
    WP1[k++] = -1.;
    WP1[k++] = -1.;

    WP1[k++] = P.x;
    WP1[k++] = 1;
    WP1[k++] = 0;

    WP1[k++] = P.y;
    WP1[k++] = 0;
    WP1[k++] = 1;
    return k;
}


//int SetWP2(const R3 &P, double *WP2) {
//    int k = 0;
//    double l[] = {1 - P.x - P.y, P.x, P.y};
//    R3 Dl[3] = {R3(-1, -1), R3(1, 0), R3(0, 1)};
//    for (int s = 0; s < 3; ++s) {
//        WP2[k++] = l[s] * (2 * l[s] - 1);
//        R3 D2 = Dl[s] * ((2 * l[s] - 1) + 2 * l[s]);
//        WP2[k++] = D2.x;
//        WP2[k++] = D2.y;
//    }
//    for (int s = 0; s < 3; ++s) {
//        int i = (s + 1) % 3, j = (s + 2) % 3;
//        WP2[k++] = 4. * l[i] * l[j];
//        R3 D2 = (4 * l[i]) * Dl[j] + (4 * l[j]) * Dl[i];
//        WP2[k++] = D2.x;
//        WP2[k++] = D2.y;
//    }
//    return k;
//}

//int SetWP0(const QF3 &qf, vector<double> &WP0, vector<double> &D) {  // Element P0
//    WP0.resize(qf.n * 3);
//    D.resize(qf.n);
//    for (int p = 0, k = 0; p < qf.n; ++p) {
//        D[p] = qf[p].a;
//        WP0[k++] = 1;
//        WP0[k++] = 0;
//        WP0[k++] = 0;
//    }
//    return 1;
//}

void SetWP1(const QF3 &qf, vector<double>& WP1, vector<double>& D) {

    for (int p = 0; p < qf.n; p++) {

        D.push_back(qf[p].a);
        R3 P = qf[p];

        WP1.push_back(1 - P.x - P.y - P.z);
        WP1.push_back(-1);
        WP1.push_back(-1);
        WP1.push_back(-1);

        WP1.push_back(P.x);
        WP1.push_back(1);
        WP1.push_back(0);
        WP1.push_back(0);

        WP1.push_back(P.y);
        WP1.push_back(0);
        WP1.push_back(1);
        WP1.push_back(0);

        WP1.push_back(P.z);
        WP1.push_back(0);
        WP1.push_back(0);
        WP1.push_back(1);
    }
}

//int SetWP2(const QF3 &qf, vector<double> &WP2, vector<double> &D) {   // Element P2
//    WP2.resize(qf.n * 6 * 3);
//    D.resize(qf.n);
//    for (int p = 0, k = 0; p < qf.n; ++p) {
//        R3 P = qf[p];
//        D[p] = qf[p].a;
//        double l[] = {1 - P.x - P.y, P.x, P.y};
//        R3 Dl[3] = {R3(-1, -1), R3(1, 0), R3(0, 1)};
//        for (int s = 0; s < 3; ++s) {
//            WP2[k++] = l[s] * (2 * l[s] - 1);
//            R3 D2 = Dl[s] * ((2 * l[s] - 1) + 2 * l[s]);
//            WP2[k++] = D2.x;
//            WP2[k++] = D2.y;
//        }
//        for (int s = 0; s < 3; ++s) {
//            int i = (s + 1) % 3, j = (s + 2) % 3;
//            WP2[k++] = 4. * l[i] * l[j];
//            R3 D2 = (4 * l[i]) * Dl[j] + (4 * l[j]) * Dl[i];
//            WP2[k++] = D2.x;
//            WP2[k++] = D2.y;
//        }
//    }
//    return 6;
//}

void SetWK(const Simplex3& simplex, const vector<double> &WP, vector<double> &WPK) {

    R3 g[4];
    simplex.GradLambdaK(g);
    int n = WP.size();

    for (int k = 0; k < n; k += simplex.vertexCount) {
        WPK[k] = WP[k];
        WPK[k + 1] = g[1].x * WP[k + 1] + g[2].x * WP[k + 2] + g[3].x * WP[k + 3];
        WPK[k + 2] = g[1].y * WP[k + 1] + g[2].y * WP[k + 2] + g[3].y * WP[k + 3];
        WPK[k + 3] = g[1].z * WP[k + 1] + g[2].z * WP[k + 2] + g[3].z * WP[k + 3];
    }
}

void SetDK(const Simplex3 &K, const vector<double> &D, vector<double> &DK) {
    // Passez de element de reference a K
    myscal(D.size(), K.volume, myscopy(D.size(), &D[0], &DK[0]));
}

void AddMass(vector<ItemOP> &vop, double alpha, int ccu, int ccv, double *vv) {
    if (ccv < 0) ccv = ccu;
    ItemOP op = {alpha, 0, 0, ccu, ccv, vv};
    vop.push_back(op);
    if (debug) cout << " AddMass " << ccu << " " << ccv << " ::: " << op << endl;

}

void AddLap(vector<ItemOP> &vop, double beta, int ccu, int ccv, double *vv) {
    if (ccv < 0) ccv = ccu;

    ItemOP op1 = {beta, 1, 1, ccu, ccv, vv};
    ItemOP op2 = {beta, 2, 2, ccu, ccv, vv};
    vop.push_back(op1);
    vop.push_back(op2);
    if (debug) cout << " AddLap " << ccu << " " << ccv << " ::: " << op1 << " " << op2 << endl;
}

void AddPdivV(vector<ItemOP> &vop, double beta, int cp, int cv1, int cv2, double *v0, double *v1) {
    ItemOP op1 = {beta, 0, 1, cp, cv1, v0};
    ItemOP op2 = {beta, 0, 2, cp, cv2, v1};
    vop.push_back(op1);
    vop.push_back(op2);
    if (debug) cout << " AddPdivV " << cp << ";" << cv1 << "," << cv2 << " ::: " << op1 << " " << op2 << endl;
}

void AdddivUQ(vector<ItemOP> &vop, double beta, int cu1, int cu2, int cq, double *v0, double *v1) {
    ItemOP op1 = {beta, 1, 0, cu1, cq, v0};
    ItemOP op2 = {beta, 2, 0, cu2, cq, v0};
    vop.push_back(op1);
    vop.push_back(op2);
    if (debug) cout << " AdddivUQ " << cu1 << "," << cu2 << ";" << cq << " ::: " << op1 << " " << op2 << endl;
}

void AddMatSymConst(hashMatrix &A, int npq, int ndofKi, const int *iK,
                    const vector<ItemOP> &vop, const double *D,
                    const double *Wi, int ndofKj, const int *jK, const double *Wj, int oi, int oj) {
    if (Wj == 0) Wj = Wi;
    if (jK == 0) jK = iK;
    if (ndofKj == 0) ndofKj = ndofKi;
    int ndofKi3 = ndofKi * 3;
    int ndofKj3 = ndofKj * 3;
    for (int ip = 0; ip < ndofKi; ++ip)
        for (int jp = 0; jp < ndofKj; ++jp) {
            double aKij = 0;
            for (int op = 0; op < vop.size(); ++op) {
                int ii = ip * 3 + vop[op].opv;
                int jj = jp * 3 + vop[op].opu;
                double cp = vop[op].c;
                double *pc = vop[op].v;

                for (int p = 0, ki = ii, kj = jj; p < npq; ++p, ki += ndofKi3, kj += ndofKj3) {
                    if (pc) {
                        cp = *pc++;
                        if (pdebugF && (ip == 0) && (jp == 0))
                            cout << p << " F[ " << pc - pdebugF - 1 << "] " << cp << endl;
                    }
                    aKij += D[p] * cp * Wi[ki] * Wj[kj];

                }
            }
            int i = iK[ip] + oi;
            int j = jK[jp] + oj;
            if (debug > 2)
                cout << "Akij " << i << " " << j << " " << aKij << endl;
            A(i, j) += aKij;
        }

}

double Integral(const Mesh3d &Th, const vector<double> &D, const vector<double> &F) {
    double r = 0;
    vector<double> DK(D);
    assert(F.size() == D.size() * Th.elementCount);
    for (int k = 0, kf = 0; k < Th.elementCount; ++k) {
        SetDK(Th[k], D, DK);
        for (int p = 0; p < DK.size(); ++p)
            r += DK[p] * F[kf++];
    }
    return r;
}

void Add2MatLap2QF(const vector<ItemOP> &vop, hashMatrix &A, const vector<double> **pWh, const vector<double> Dh,
                   const vector<int> **pnum, int *offset, const Mesh3d &Th) {
    long npq = Dh.size();
    typedef map<pair<int, int>, vector<ItemOP> > MOPB;
    typedef MOPB::const_iterator CI_MOPB;
    MOPB opb;

    for (int i = 0; i < vop.size(); ++i)
        opb[make_pair(vop[i].cu, vop[i].cv)].push_back(vop[i]);

    for (CI_MOPB ib = opb.begin(); ib != opb.end(); ++ib) {
        //  block
        int cu = ib->first.first;
        int cv = ib->first.second;
        vector<ItemOP> vopuv = ib->second;
        // to store the adresse of pointeur of function quadrature
        vector<double **> pvopuv;
        for (int i = 0; i < vopuv.size(); ++i)
            if (vopuv[i].v) pvopuv.push_back(&(vopuv[i].v)); // on stoke les adr pointeur
        long nptr = pvopuv.size();
        const vector<int> &numi = *pnum[cv];
        const vector<int> &numj = *pnum[cu];
        const vector<double> &Whi = *pWh[cv];
        const vector<double> &Whj = *pWh[cu];

        int oi = offset[cv];
        int oj = offset[cu];
        long ndofKi = numi.size() / Th.elementCount;
        long ndofKj = numj.size() / Th.elementCount;

        vector<double> WKi(Whi.size()), WKj(Whj.size()), DK(Dh.size());
        cout << "\n ** " << cu << " " << cv << "  ; " << oi << " " << oj << " || " << ndofKi << " " << ndofKj
             << " :  " << vopuv << endl;

        for (int k = 0; k < Th.elementCount; ++k) {
            const Simplex3 &K = Th[k];
            SetWK(K, Whi, WKi);
            SetWK(K, Whj, WKj);
            SetDK(K, Dh, DK);
            AddMatSymConst(A, npq, ndofKi, &numi[k * ndofKi], vopuv, &DK[0], &WKi[0],
                           ndofKj, &numj[k * ndofKj], &WKj[0], oi, oj
            );
            for (int ptr = 0; ptr < nptr; ++ptr) {
                *(pvopuv[ptr]) += npq; // decalage de pointeur dans dans le item ... en fin de boucle
                if (pdebugF) cout << ptr << " " << k << " " << *(pvopuv[ptr]) - pdebugF << " " << npq << endl;
            }

        }
        //      for(int ptr=0; ptr < nptr; ++ptr)
        //          *(pvopuv[ptr]) += npq; // decalage de pointeur dans dans les items en fin de boucle  ...

    }
}

void Add2MatLap2QF(hashMatrix &A, const vector<double> &Wh, const vector<double> &Dh, vector<ItemOP> vop,
                   const vector<int> &num, const Mesh3d &Th) {
    int npq = Dh.size();
    int ndofK = num.size() / Th.elementCount;

    vector<double **> pvop;

    for (int i = 0; i < vop.size(); ++i)
        if (vop[i].v) {
            pvop.push_back(&(vop[i].v));
        }// on stocke les adr pointeur des fonctions quadrature
    int nptr = pvop.size();

    vector<double> WK(Wh.size()), DK(Dh.size());

    for (int k = 0; k < Th.elementCount; ++k) {
        const Simplex3 &K = Th[k];
        SetWK(K, Wh, WK);
        SetDK(K, Dh, DK);
        AddMatSymConst(A, npq, ndofK, &num[k * ndofK], vop, &DK[0], &WK[0]);
        for (int ptr = 0; ptr < nptr; ++ptr) {
            *(pvop[ptr]) += npq; // decalage de pointeur dans dans le item ... en fin de boucle
            if (pdebugF) cout << ptr << " " << k << " " << *(pvop[ptr]) - pdebugF << " " << npq << endl;
        }

    }

}

double min(double *v, int n) {
    double z0 = v[0];
    for (int i = 1; i < n; ++i)
        z0 = min(z0, v[i]);
    return z0;
}

double max(double *v, int n) {
    double z0 = v[0];
    for (int i = 1; i < n; ++i)
        z0 = max(z0, v[i]);
    return z0;
}

double normsum(double *u, double *v, int n) {
    double err = 0;
    for (int i = 0; i < n; ++i)
        err = max(err, abs(v[i] - u[i]));
    return err;
}

int buildIndexes(const Mesh3d &Th, int ndofs, int ndofa, int ndoft, vector<int> &num) {
    int ndofK = ndofs * 3 + ndofa * 3 + ndoft;
    num.resize(ndofK * Th.elementCount);
    int dof = 0;
    int knum = 0;
    map<pair<int, int>, int> M;
    for (int k = 0; k < Th.elementCount; ++k) {
        if (ndofs)
            for (int sp = 0; sp < 3; ++sp) {
                int s = Th(k, sp);
                pair<int, int> key(s, -1);
                if (M.find(key) == M.end()) {
                    M[key] = dof;
                    dof += ndofs;
                }
                int sdof = M[key];
                for (int i = 0; i < ndofs; ++i)
                    num[knum++] = sdof + i;

            }
        if (ndofa)
            for (int a = 0; a < 3; ++a) {// Arete oppose au sommet
                int s1 = Th(k, (a + 1) % 3);
                int s2 = Th(k, (a + 2) % 3);
                if (s2 > s1) std::swap(s1, s2);
                pair<int, int> key(s1, s2);
                if (M.find(key) == M.end()) {
                    M[key] = dof;
                    dof += ndofa;
                }
                int adof = M[key];
                for (int i = 0; i < ndofa; ++i)
                    num[knum++] = adof + i;

            }
        if (ndoft) {//
            pair<int, int> key(k, -2);
            if (M.find(key) == M.end()) {
                M[key] = dof;
                dof += ndoft;
            }
            int tdof = M[key];
            for (int i = 0; i < ndoft; ++i)
                num[knum++] = tdof + i;
        }

    }
    assert(knum == num.size());
    return dof;
}

ostream &operator<<(ostream &f, const ItemOP &o) {
    static string so[3] = {" ", ".dx ", ".dy "};
    return f << std::showpos << " " << o.c << std::noshowpos << " u" << o.cu + 1 << so[o.opu]
             << "v" << o.cv + 1 << so[o.opv] << " ";
}

ostream &operator<<(ostream &f, const vector<ItemOP> &vop) {
    for (int i = 0; i < vop.size(); ++i)
        f << vop[i];
    return f;
}

//void SetF3(const Mesh3d &Th, const QF3 &qf, const vector<double> &WP2,
//           const vector<int> &num, const vector<double> &u, vector<double> &F3,
//           void (*FF)(double *f, double x, double y, double v, double v_x, double v_y)
//) {
//
//    vector<double> WK(WP2.size());
//    int lF3 = F3.size();
//    int ndofK = num.size() / Th.elementCount;
//    assert(lF3 == ndofK * 3 * qf.n);
//
//    for (int k = 0, kf = 0, nk = 0; k < Th.elementCount; ++k, nk += ndofK) {
//        const Simplex3 &K = Th[k];
//        SetWK(K, WP2, WK);
//
//        for (int p = 0, ii = 0; p < qf.n; ++p, kf += 3) {
//            R3 P = K(qf[p]);
//            double uv = 0, ux = 0, uy = 0;
//            for (int ip = 0; ip < ndofK; ++ip, ii += 3) {
//                uv += u[num[nk + ip]] * WK[ii];
//                ux += u[num[nk + ip]] * WK[ii + 1];
//                uy += u[num[nk + ip]] * WK[ii + 2];
//            }
//            FF(&F3[kf], P.x, P.y, uv, ux, uy);
//            if (debug || pdebugF)
//                cout << k << " : F " << p << " kf " << kf << " : "
//                     << F3[kf] << " " << F3[kf + 1] << " " << F3[kf + 2] << endl;
//        }
//    }
//}

//void SetF(const Mesh3d &Th, const QF3 &qf, const vector<double> &WP2, const vector<int> &num, const vector<double> &u,
//          vector<double> &F, double (*FF)(double x, double y, double v, double v_x, double v_y)) {
//
//    vector<double> WK(WP2.size());
//    int ndofK = num.size() / Th.elementCount;
//    // int nwk = ndofK*3*qf.n;
//    for (int k = 0, kf = 0, nk = 0; k < Th.elementCount; ++k, nk += ndofK) {
//        const Simplex3 &K = Th[k];
//        SetWK(K, WP2, WK);
//
//        for (int p = 0, ii = 0; p < qf.n; ++p, ++kf) {
//            R3 P = K(qf[p]);
//            double uv = 0, ux = 0, uy = 0;
//            for (int ip = 0; ip < ndofK; ++ip, ii += 3) {
//                uv += u[num[nk + ip]] * WK[ii];
//                ux += u[num[nk + ip]] * WK[ii + 1];
//                uy += u[num[nk + ip]] * WK[ii + 2];
//            }
//            F[kf] = FF(P.x, P.y, uv, ux, uy);
//            if (debug || pdebugF)
//                cout << k << " : F " << p << " kf " << kf << "  == " << F[kf] << endl;
//        }
//    }
//}
