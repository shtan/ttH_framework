#ifndef commonfunc_h
#define commonfunc_h

#include <vector>
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "Math/QuantFuncMathCore.h"

using namespace std;

namespace commonfunc
{

//Functions to convert between polar and Cartesian coords
inline TLorentzVector lorentz_maker_pol(double pt, double eta, double phi, double m) {
    TLorentzVector vec;
    vec.SetPtEtaPhiM(pt, eta, phi, m);
    return vec;
}

inline TLorentzVector lorentz_maker_car(double px, double py, double pz, double e){
    TLorentzVector vec;
    vec.SetPxPyPzE(px, py, pz, e);
    return vec;
}

inline TLorentzVector lorentz_maker_carM(double px, double py, double pz, double m){
    double e = sqrt( pow(px, 2) + pow(py, 2) + pow(pz, 2) + pow(m, 2) );
    TLorentzVector vec;
    vec.SetPxPyPzE(px, py, pz, e);
    return vec;
}

//Define breitWigner function to calculate mTop and mW chi2
inline double breitWignerErr(const double &mass, const double &width, const double &deltaMass)
{
    double scaledDeltaMass = deltaMass * width;
    double scaledDeltaMass2 = scaledDeltaMass * scaledDeltaMass;
    double Zscore = ROOT::Math::normal_quantile(
                0.31830988618379067154 *atan2(scaledDeltaMass2 + 2. * scaledDeltaMass * mass,
                mass * width) + 0.5, 1.0);
    return Zscore * Zscore;
}

}

#endif
