#include "WDaughterEllipseCalculator.h"
#include <complex>

using namespace commonstruct;

//WDaughterEllipseCalculator::WDaughterEllipseCalculator(top_sys& topsystem) {}
WDaughterEllipseCalculator::WDaughterEllipseCalculator( top_system& topsystem, int &debug )
    : topsys(topsystem),
      Ab_(4, 4), AWDaughter1_(4, 4), Htilde_(3, 3), H_(3, 3), Hperp_(3, 3),
      HperpInv_(3, 3), Nperp_(3, 3), WDaughterPerp_(3), pWDaughter_(3),
      errorFlag_(false), debug_verbosity(debug)
{
    if (debug_verbosity >= 2)
        cout << "Starting WDaughterEllipseCalculator constructor" << endl;
    
}
 


WDaughterEllipseCalculator::~WDaughterEllipseCalculator()
{
    // cout << "destructor" << endl;
}

void WDaughterEllipseCalculator::preSetupEllipse()
{
    setBJetFactors();
    setMeasuredWDaughterFactors();
    setAngles();
    initializeMatrices();
}

void WDaughterEllipseCalculator::setupEllipsePart2()
{
    Wsurface();
    measuredWDaughterEllipsoid();
    bJetEllipsoid();
    calcZ2();
}


double WDaughterEllipseCalculator::getZ2()
{
    setupEllipsePart2();
    // cout << "Top mass is " <<topsys.calc.mTop() << endl;
    // cout << "W mass is " << topsys.calc.mW() << endl;
    // cout << "second W daughter mass is " << topsys.calc.Wd2_m() << endl;
    return Z2_;
}


void WDaughterEllipseCalculator::setBJetFactors()
{
    // cout << "setting b-jet factors" << endl;
    bJetP2_ = (topsys.calc.b_px() * topsys.calc.b_px() + topsys.calc.b_py() * topsys.calc.b_py() + topsys.calc.b_pz() * topsys.calc.b_pz());
    bJetP_ = sqrt(bJetP2_);
    double bJetE2 = topsys.calc.b_e() * topsys.calc.b_e();
    //bJetMass2_ = bJetE2 - bJetP2_;
    bJetMass2_ = topsys.calc.b_m() * topsys.calc.b_m();
    bJetBeta2_ = bJetP2_ / bJetE2;
    bJetBeta_ = sqrt(bJetBeta2_);
    bJetGamma2_ = 1.0 / (1.0 - bJetBeta2_);
    bJetGamma_ = sqrt(bJetGamma2_);
}

void WDaughterEllipseCalculator::setMeasuredWDaughterFactors()
{
    // cout << "setting measured W daughter factors" << endl;
    WDaughter1P2_ =
        (topsys.calc.Wd1_px() * topsys.calc.Wd1_px() + topsys.calc.Wd1_py() * topsys.calc.Wd1_py() +
         topsys.calc.Wd1_pz() * topsys.calc.Wd1_pz());
    WDaughter1P_ = sqrt(WDaughter1P2_);
    double WDaughter1E2 = topsys.calc.Wd1_e() * topsys.calc.Wd1_e();
    //WDaughter1Mass2_ = max(0., WDaughter1E2 - WDaughter1P2_);
    WDaughter1Mass2_ = topsys.calc.Wd1_m() * topsys.calc.Wd1_m();
    WDaughter1Beta2_ = WDaughter1P2_ / WDaughter1E2;
    WDaughter1Beta_ = sqrt(WDaughter1Beta2_);
    WDaughter1Gamma2_ = 1.0 / (1.0 - WDaughter1Beta2_);
    WDaughter1Gamma_ = sqrt(WDaughter1Gamma2_);
    //WDaughter1Phi_ = (topsys.calc.Wd1_px() == 0.0 && topsys.calc.Wd1_py() == 0.0)
    //                     ? 0.0
    //                     : atan2(topsys.calc.Wd1_py(), topsys.calc.Wd1_px());
    WDaughter1Theta_ =
        (topsys.calc.Wd1_px() == 0.0 && topsys.calc.Wd1_py() == 0.0 && topsys.calc.Wd1_pz() == 0.0)
            ? 0.0
            : atan2(sqrt(topsys.calc.Wd1_px() * topsys.calc.Wd1_px() +
                         topsys.calc.Wd1_py() * topsys.calc.Wd1_py()),
                    topsys.calc.Wd1_pz());
}

void WDaughterEllipseCalculator::setAngles()
{
    // cout << "setting angles" << endl;
    double WDaughter1DotbJet =
        (topsys.calc.Wd1_px() * topsys.calc.b_px() + topsys.calc.Wd1_py() * topsys.calc.b_py() +
         topsys.calc.Wd1_pz() * topsys.calc.b_pz());
    c2_ = WDaughter1DotbJet * WDaughter1DotbJet / (bJetP2_ * WDaughter1P2_);
    c_ = WDaughter1DotbJet / sqrt(bJetP2_ * WDaughter1P2_);
    s2_ = 1. - c2_;
    s_ = sqrt(s2_);
    // cout << "Cosine: " << c_ << endl;
    // cout << "Sine  : " << s_ << endl;
}

void WDaughterEllipseCalculator::initializeMatrices()
{
    // cout << "initializing matrices" << endl;

    Ab_.Zero();
    AWDaughter1_.Zero();
    Htilde_.Zero();
    H_.Zero();
    Hperp_.Zero();
    HperpInv_.Zero();
    Nperp_.Zero();

    WDaughterPerp_.Zero();
    pWDaughter_.Zero();

    // cout << "the b-jet matrix has " << Ab_.GetNrows() << " rows and " <<
    // Ab_.GetNcols() << " columns"  << endl;
    // cout << "the measuredWDaughter matrix has " << AWDaughter1_.GetNrows() <<
    // " rows and " << AWDaughter1_.GetNcols() << " columns"  << endl;
}

void WDaughterEllipseCalculator::Wsurface()
{
    // cout << "calculating the W surface" << endl;
    x0p_ = -(0.5 / topsys.calc.b_e()) * (topsys.calc.mTop2() - topsys.calc.mW2() - bJetMass2_);
    x0_ = -(0.5 / topsys.calc.Wd1_e()) * (topsys.calc.mW2() - topsys.calc.Wd2_m2() - WDaughter1Mass2_);
    Sx_ =
        x0_ / WDaughter1Beta_ - WDaughter1P_ / WDaughter1Beta2_ + WDaughter1P_;
    epsilon2_ = (topsys.calc.mW2() - topsys.calc.Wd2_m2()) -
                WDaughter1Beta2_ * (topsys.calc.mW2() - topsys.calc.Wd2_m2());
    // cout << "b-jet energy is " << bJetE_ << endl;
    // cout << "b-jet mass^2 is " << bJetMass2_ << endl;
    // cout << "mW^2 is " << topsys.calc.mW2() << endl;
    // cout << "mTop^2 is " << topsys.calc.mTop2() << endl;
    // cout << "x0p is " << x0p_ << endl;
    // cout << "1st W daughter energy is " << topsys.calc.Wd1_e() << endl;
    // cout << "1st W daughter mass^2 is " << WDaughter1Mass2_ << endl;
    // cout << "2nd W daughter mass^2 is " << topsys.calc.Wd2_m2() << endl;
    // cout << "x0 is " << x0_ << endl;
    // cout << "c is " << c_ << endl;
    // cout << "s is " << s_ << endl;
    // cout << "Sx is " << Sx_ << endl;
    // cout << "epsilon^2 is " << epsilon2_ << endl;
}

void WDaughterEllipseCalculator::bJetEllipsoid()
{
     //cout << "calculating the b-jet ellipsoid" << endl;

    Ab_[0][0] = 1 - c2_ * bJetBeta2_;
    Ab_[1][0] = -c_ * s_ * bJetBeta2_;
    Ab_[2][0] = 0;
    Ab_[3][0] = c_ * x0p_ * bJetBeta_;

    Ab_[0][1] = -c_ * s_ * bJetBeta2_;
    Ab_[1][1] = 1 - s2_ * bJetBeta2_;
    Ab_[2][1] = 0;
    Ab_[3][1] = s_ * x0p_ * bJetBeta_;

    Ab_[0][2] = 0;
    Ab_[1][2] = 0;
    Ab_[2][2] = 1;
    Ab_[3][2] = 0;

    Ab_[0][3] = c_ * x0p_ * bJetBeta_;
    Ab_[1][3] = s_ * x0p_ * bJetBeta_;
    Ab_[2][3] = 0;
    Ab_[3][3] = topsys.calc.mW2() - x0p_ * x0p_;

    // cout << "Measured b-jet ellipsoid:" << endl;
    // Ab_.Print();
}

void WDaughterEllipseCalculator::measuredWDaughterEllipsoid()
{
    // cout << "calculating the measured W daughter ellipsoid" << endl;

    AWDaughter1_[0][0] = 1. - WDaughter1Beta2_;
    AWDaughter1_[1][0] = 0;
    AWDaughter1_[2][0] = 0;
    AWDaughter1_[3][0] = Sx_ * WDaughter1Beta2_;

    AWDaughter1_[0][1] = 0;
    AWDaughter1_[1][1] = 1;
    AWDaughter1_[2][1] = 0;
    AWDaughter1_[3][1] = 0;

    AWDaughter1_[0][2] = 0;
    AWDaughter1_[1][2] = 0;
    AWDaughter1_[2][2] = 1;
    AWDaughter1_[3][2] = 0;

    AWDaughter1_[0][3] = Sx_ * WDaughter1Beta2_;
    AWDaughter1_[1][3] = 0;
    AWDaughter1_[2][3] = 0;
    AWDaughter1_[3][3] = topsys.calc.mW2() - x0_ * x0_ - epsilon2_;

    // cout << "Measured W daughter ellipsoid:" << endl;
    // AWDaughter1_.Print();
}

void WDaughterEllipseCalculator::calcZ2()
{
    // cout << "Calculating Z^2" << endl;
    Sy_ = (1. / s_) * (x0p_ / bJetBeta_ - c_ * Sx_);
    // cout << "Sy is " << Sy_ << endl;
    omega_ = (1. / s_) *
             (WDaughter1Beta_ / bJetBeta_ - c_); // only the positive slope
    // cout << "omega is " << omega_ << endl;
    double Omega2 = omega_ * omega_ + 1. - WDaughter1Beta2_;
    // cout << "Omega^2 is " << Omega2 << endl;
    Omega_ = sqrt(Omega2);
    // cout << "Omega is " << Omega_ << endl;
    x1_ = Sx_ - (1. / Omega2) * (Sx_ + omega_ * Sy_);
    y1_ = Sy_ - (1. / Omega2) * omega_ * (Sx_ + omega_ * Sy_);
    Z2_ = x1_ * x1_ * Omega2 - (Sy_ - omega_ * Sx_) * (Sy_ - omega_ * Sx_) -
          (topsys.calc.mW2() - x0_ * x0_ - epsilon2_);
    // if(Z2_ > 0) cout << "Z^2 is positive: " << Z2_ << endl;
    // else if(Z2_ < 0) cout << "Z^2 is negative: " << Z2_ << endl;
    // else if(Z2_ == 0) cout << "Z^2 is exactly zero" << endl;
    // else cout << "Z^2 is nan" << endl;
    // if(Z2_<0.) Htilde_.Print();
}

void WDaughterEllipseCalculator::WDaughterSolution()
{
    // cout << "calculating second W daughter ellipse" << endl;
    // cout << "Current top mass is " << topsys.calc.mTop() << endl;

    calcZ2();

    // cout << "checking sign of Z^2: " << Z2_ << endl;

    if (Z2_ < 0) {
        // cout << "Z^2 is negative !!" << endl;
        // cout << "Z^2 = " << Z2_ << endl;

        // Vary the top mass to force Z^2>0
        // x0p_ depends on topsys.calc.mTop(), so Wsurface and bJetEllipsoid need to be
        // called to propagate the changes to the top mass

        errorFlag_ = true;
        topsys.vars.error_flag = true;
        Z2_ = 0;

        //      calcTopMassCorrection(); //calls to Wsurface and bJetEllipsoid
        //
        //      if( Z2_<0 ) Z2_=0; //for now reject events

    } else {
        // cout << "Z^2 is positive, no need to vary the top mass" << endl;
        errorFlag_ = false;
        topsys.vars.error_flag = false;
    }

    //topsys.vars.error_flag = errorFlag_;

    double Z = sqrt(Z2_);

    Htilde_[0][0] = Z / Omega_;
    Htilde_[0][1] = 0;
    Htilde_[0][2] = x1_ - WDaughter1P_;

    Htilde_[1][0] = omega_ * Z / Omega_;
    Htilde_[1][1] = 0;
    Htilde_[1][2] = y1_;

    Htilde_[2][0] = 0;
    Htilde_[2][1] = Z;
    Htilde_[2][2] = 0;

    // Htilde_.Print();
}

TMatrixD WDaughterEllipseCalculator::rotationMatrix(int axis, double angle)
{
    TMatrixD r(3, 3);
    r.Zero();
    if (axis != 0 && axis != 1 && axis != 2)
        return r;

    for (int i = 0; i < 3; i++) {
        r[i][i] = cos(angle);
    }

    for (int i = -1; i <= 1; i++) {
        double row = (axis - i) % 3;
        if (row < 0)
            row += 3;
        double col = (axis + i) % 3;
        if (col < 0)
            col += 3;
        r[row][col] = i * sin(angle) + (1 - i * i);
    }

    return r;
}

void WDaughterEllipseCalculator::labSystemTransform()
{
    // cout << "Boosting back into the lab frame" << endl;

    // rotate Htilde to H
    TMatrixD Rz = rotationMatrix(2, - topsys.calc.Wd1_phi());
    TMatrixD Ry = rotationMatrix(1, 0.5 * M_PI - WDaughter1Theta_);
    double bJetP[3] = {topsys.calc.b_px(), topsys.calc.b_py(), topsys.calc.b_pz()};
    TMatrixD bJet_xyz(3, 1, bJetP);
    TMatrixD rM(Ry, TMatrixD::kMult, TMatrixD(Rz, TMatrixD::kMult, bJet_xyz));
    double *rA = rM.GetMatrixArray();
    double phi = -TMath::ATan2(rA[2], rA[1]);
    TMatrixD Rx = rotationMatrix(0, phi);

    H_ = TMatrixD(TMatrixD::kTransposed, Rz);
    TMatrixD RyT(TMatrixD::kTransposed, Ry);
    TMatrixD RxT(TMatrixD::kTransposed, Rx);

    H_ *= RyT;
    H_ *= RxT;
    H_ *= Htilde_;

    // calculate Hperp
    double Hvalues[9] = {H_[0][0], H_[0][1], H_[0][2], H_[1][0], H_[1][1],
                         H_[1][2], 0,        0,        1};
    TArrayD Harray(9, Hvalues);
    Hperp_.SetMatrixArray(Harray.GetArray());

    // Hperp_.Print();
}

void WDaughterEllipseCalculator::calcWDaughterEllipse()
{
    Wsurface();
    // printFactors();
    measuredWDaughterEllipsoid();
    bJetEllipsoid();
    // printFactors();
    WDaughterSolution();
    labSystemTransform();
    // Ab_.Print();
    // AWDaughter1_.Print();
    // printFactors();
}

void WDaughterEllipseCalculator::calcExtendedWDaughterEllipse()
{
    // avoid non-invertible matrices
    if (Z2_ <= 0) {
        // cout << "Z^2 is negative: " << Z2_ << ", cannot calculate extended
        // representation" << endl;
        Nperp_.Zero();
        return;
    }
    // cout << "Z^2 is positive: " << Z2_ << ", now calculating extended
    // representation" << endl;
    // cout << "Z^2 = " << Z2_ << endl;
    // Hperp_.Print();
    HperpInv_ = Hperp_;
    HperpInv_.Invert();
    TMatrixD U(3, 3);
    U.Zero();
    U[0][0] = 1;
    U[1][1] = 1;
    U[2][2] = -1;
    Nperp_ = TMatrixD(HperpInv_, TMatrixD::kTransposeMult,
                      TMatrixD(U, TMatrixD::kMult, HperpInv_));
}

TMatrixD *WDaughterEllipseCalculator::getExtendedWDaughterEllipse()
{
    return &Nperp_;
}

TMatrixD *WDaughterEllipseCalculator::getHomogeneousWDaughterEllipse()
{
    return &Hperp_;
}

void WDaughterEllipseCalculator::getWDaughterMomentum()
{
    // calcExtendedWDaughterEllipse();
    // Hperp_.Print();
    double tArray[3] = {cos(topsys.vars.theta), sin(topsys.vars.theta), 1.};
    TVectorD t(3, tArray);
    WDaughterPerp_ = t;
    WDaughterPerp_ *= Hperp_;
    pWDaughter_ = WDaughterPerp_;
    pWDaughter_ *= HperpInv_;
    pWDaughter_ *= H_;
    topsys.vars.Wd2_px = pWDaughter_[0];
    topsys.vars.Wd2_py = pWDaughter_[1];
    topsys.vars.Wd2_pz = pWDaughter_[2];
    //return &pWDaughter_;
}


