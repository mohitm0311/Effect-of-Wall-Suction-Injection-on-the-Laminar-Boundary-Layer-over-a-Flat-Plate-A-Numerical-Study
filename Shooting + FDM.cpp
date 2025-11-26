#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
using namespace std;


// Forward Euler FDM step
void FD_step(double &f, double &g, double &h, double hstep) {
    double f_new = f + hstep * g;            // f' = g
    double g_new = g + hstep * h;            // g' = h
    double h_new = h - 0.5 * hstep * f * h;  // h' = -0.5 f h

    f = f_new;
    g = g_new;
    h = h_new;
}

// Integrate once for a given guess a = f''(0)
double integrate_once_FDM(double S, double a,
                          double etamax, double hstep,
                          vector<double>* ETA = nullptr,
                          vector<double>* F = nullptr,
                          vector<double>* FP = nullptr)
{
    double f = S;    // f(0) = S
    double g = 0.0;  // f'(0) = 0
    double h = a;    // f''(0) = a  (shooting parameter)

    double eta = 0.0;

    if (ETA && F && FP) {
        ETA->clear(); F->clear(); FP->clear();
        ETA->push_back(eta);
        F->push_back(f);
        FP->push_back(g);
    }

    while (eta < etamax - 1e-12) {
        FD_step(f, g, h, hstep);
        eta += hstep;

        if (ETA && F && FP) {
            ETA->push_back(eta);
            F->push_back(f);
            FP->push_back(g);
        }
    }

    // Boundary condition: f'(∞) = g(etamax) = 1
    return g - 1.0;
}


// Secant shooting method
double shoot_FDM(double S,
                 double etamax, double hstep,
                 double tol, int maxIter,
                 vector<double> &ETA,
                 vector<double> &F,
                 vector<double> &FP)
{
    double a1 = 0.3;
    double a2 = 0.4;

    double R1 = integrate_once_FDM(S, a1, etamax, hstep);
    double R2 = integrate_once_FDM(S, a2, etamax, hstep);

    for (int it = 0; it < maxIter; it++) {
        double denom = R2 - R1;
        if (fabs(denom) < 1e-14) break;

        double a3 = a2 - R2 * (a2 - a1) / denom;
        double R3 = integrate_once_FDM(S, a3, etamax, hstep);

        if (fabs(R3) < tol) {
            integrate_once_FDM(S, a3, etamax, hstep, &ETA, &F, &FP);
            cout << "FDM converged in " << it+1 << " iterations.\n";
            return a3;
        }

        a1 = a2;  R1 = R2;
        a2 = a3;  R2 = R3;
    }

    integrate_once_FDM(S, a2, etamax, hstep, &ETA, &F, &FP);
    return a2;
}

int main() {
    double S = 0.2;         // ONLY S = 0.2 (for validation)
    double etamax = 8.0;
    double hstep  = 0.01;  
    double tol    = 1e-6;
    int maxIter   = 60;

    vector<double> ETA, F, FP;

    // Solve using Shooting + FDM
    double fpp0 = shoot_FDM(S, etamax, hstep, tol, maxIter, ETA, F, FP);


    cout << "  RESULTS FOR FDM + SHOOTING\n";
    cout << "S = 0.2\n";
    cout << "f''(0) = " << fpp0 << endl;


    ofstream out("profile_FDM_S_0.20.txt");
    out << "# eta   f(eta)   f'(eta)\n";
    for (int i = 0; i < ETA.size(); i++)
        out << ETA[i] << " " << F[i] << " " << FP[i] << "\n";
    out.close();
    return 0;
}