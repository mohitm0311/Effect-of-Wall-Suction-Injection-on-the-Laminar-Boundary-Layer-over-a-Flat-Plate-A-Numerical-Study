#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
using namespace std;

// Blasius system:
// f' = g
// g' = h
// h' = -0.5 * f * h

double fxn_f(double eta, double f, double g, double h){
    return g;
}

double fxn_g(double eta, double f, double g, double h){
    return h;
}

double fxn_h(double eta, double f, double g, double h){
    return -0.5 * f * h;
}


void RK4_step(double &f, double &g, double &h, double eta, double hstep){

    double k1_f = fxn_f(eta,       f,           g,           h);
    double k1_g = fxn_g(eta,       f,           g,           h);
    double k1_h = fxn_h(eta,       f,           g,           h);

    double f2 = f + 0.5*hstep*k1_f;
    double g2 = g + 0.5*hstep*k1_g;
    double h2 = h + 0.5*hstep*k1_h;

    double k2_f = fxn_f(eta+0.5*hstep, f2,      g2,          h2);
    double k2_g = fxn_g(eta+0.5*hstep, f2,      g2,          h2);
    double k2_h = fxn_h(eta+0.5*hstep, f2,      g2,          h2);

    double f3 = f + 0.5*hstep*k2_f;
    double g3 = g + 0.5*hstep*k2_g;
    double h3 = h + 0.5*hstep*k2_h;

    double k3_f = fxn_f(eta+0.5*hstep, f3,      g3,          h3);
    double k3_g = fxn_g(eta+0.5*hstep, f3,      g3,          h3);
    double k3_h = fxn_h(eta+0.5*hstep, f3,      g3,          h3);

    double f4 = f + hstep*k3_f;
    double g4 = g + hstep*k3_g;
    double h4 = h + hstep*k3_h;

    double k4_f = fxn_f(eta+hstep,     f4,      g4,          h4);
    double k4_g = fxn_g(eta+hstep,     f4,      g4,          h4);
    double k4_h = fxn_h(eta+hstep,     f4,      g4,          h4);

    f += (hstep/6.0)*(k1_f + 2*k2_f + 2*k3_f + k4_f);
    g += (hstep/6.0)*(k1_g + 2*k2_g + 2*k3_g + k4_g);
    h += (hstep/6.0)*(k1_h + 2*k2_h + 2*k3_h + k4_h);
}

// Integrate once for a given S and guess a = f''(0)
// If E,F,G are given -> store profile
double integrate_once(double S, double a,
                      double etamax, double hstep,
                      vector<double>* E = nullptr,
                      vector<double>* F = nullptr,
                      vector<double>* G = nullptr)
{
    double f = S;     // f(0) = S  (suction/injection parameter)
    double g = 0.0;   // f'(0) = 0 (no-slip)
    double h = a;     // f''(0) = a (shooting parameter)

    double eta = 0.0;

    if(E && F && G){
        E->clear(); F->clear(); G->clear();
        E->push_back(eta);
        F->push_back(f);
        G->push_back(g);
    }

    while(eta < etamax - 1e-12){
        RK4_step(f, g, h, eta, hstep);
        eta += hstep;

        if(E && F && G){
            E->push_back(eta);
            F->push_back(f);
            G->push_back(g);
        }
    }

    return g - 1.0;
}

// Shooting method using secant 
double shoot(double S,
             double etamax, double hstep,
             double tol, int maxIter,
             vector<double> &E,
             vector<double> &F,
             vector<double> &G)
{
    // initial guesses for f''(0),
   
    double a1 = 0.3;
    double a2 = 0.4;

    double R1 = integrate_once(S, a1, etamax, hstep);
    double R2 = integrate_once(S, a2, etamax, hstep);

    for(int it=0; it<maxIter; it++){

        double denom = (R2 - R1);
        if(fabs(denom) < 1e-14) break; 

        double a3 = a2 - R2*(a2 - a1)/denom;
        double R3 = integrate_once(S, a3, etamax, hstep);

        if(fabs(R3) < tol){
            integrate_once(S, a3, etamax, hstep, &E, &F, &G);
            cout<<"S = "<<S<<" converged in "<<it+1<<" iterations\n";
            return a3;
        }

       
        a1 = a2; R1 = R2;
        a2 = a3; R2 = R3;
    }

    // if not converged, just output last profile for a2
    cout<<"S = "<<S<<" did not fully converge, |R| = "<<fabs(R2)<<endl;
    integrate_once(S, a2, etamax, hstep, &E, &F, &G);
    return a2;
}

int main(){
    double etamax = 8.0;
    double hstep  = 0.01;
    double tol    = 1e-6;
    int maxIter = 40;

    // choose S values: suction (>0), injection (<0), classical (0)
    vector<double> Svals = {-0.5, -0.2, 0.0, 0.2, 0.5};

    // Cf = 2 f''(0)/sqrt(Re_x)
    double Rex = 1.0e5; // example Reynolds Number

    for(double S : Svals){

        vector<double> eta, f, fp;
        double fpp0 = shoot(S, etamax, hstep, tol, maxIter, eta, f, fp);

        double Cf = 2.0 * fpp0 / sqrt(Rex);

        cout<<"S = "<<S<<"  f''(0) = "<<fpp0<<"  Cf = "<<Cf<<endl;

        // write file for plotting
        string name = "profile_S_" + to_string(S) + ".txt";
        ofstream out(name.c_str());
        out<<"# eta   f(eta)   f'(eta)\n";
        for(int i=0; i< (int)eta.size(); i++){
            out<<eta[i]<<" "<<f[i]<<" "<<fp[i]<<"\n";
        }
        out.close();
    }

    return 0;
}
