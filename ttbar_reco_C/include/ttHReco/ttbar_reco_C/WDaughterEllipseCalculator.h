// C++ implementation of arXiv:1305.1872
// Analytic solutions for constrained momentum of second W daughter in decay of
// top quarks
// Gala Nicolas Kaufman (gnn4@cornell.edu)

#include <memory>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include <list>
#include <utility>

#include <TVectorD.h>
#include <TMatrixD.h>
#include <TArrayD.h>
#include <TMath.h>
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

#include "Math/RootFinderAlgorithms.h"
#include "Math/Polynomial.h"

#include "commonstruct.h"

using namespace std;
using namespace ROOT::Math;
//using namespace commonstruct;

class WDaughterEllipseCalculator
{

  private:
    commonstruct::top_system &topsys;

    double bJetP2_, bJetP_, bJetMass2_;
    double bJetBeta_, bJetBeta2_, bJetGamma_, bJetGamma2_;
    double WDaughter1Beta_, WDaughter1Beta2_, WDaughter1Gamma_,
        WDaughter1Gamma2_;
    double WDaughter1P2_, WDaughter1P_, WDaughter1Mass2_;
    //double WDaughter1Phi_;
    double WDaughter1Theta_;

    // parameters
    double x0_, x0p_;
    double Sx_, Sy_;
    double epsilon2_;

    double c_, s_, c2_, s2_; // cosine and sine of theta_{b,l}

    double omega_;
    double Omega_;
    double x1_, y1_;
    double Z2_;

    // matrices
    TMatrixD Ab_;
    TMatrixD AWDaughter1_;

    TMatrixD Htilde_;
    TMatrixD H_;
    TMatrixD Hperp_;
    TMatrixD HperpInv_;

    TMatrixD Nperp_;

    TVectorD WDaughterPerp_;
    TVectorD pWDaughter_;

    int& debug_verbosity;

  public:
    void setBJetFactors();
    void setMeasuredWDaughterFactors();

    double getZ2();

    void setAngles();

    void initializeMatrices();

    TMatrixD rotationMatrix(int, double);

    void Wsurface();
    void bJetEllipsoid();
    void measuredWDaughterEllipsoid();
    void calcZ2();
    void WDaughterSolution();
    void labSystemTransform();

    bool errorFlag_;

    WDaughterEllipseCalculator( commonstruct::top_system&, int& );
    ~WDaughterEllipseCalculator();

    bool badPoint() { return errorFlag_; };

    void preSetupEllipse();
    void setupEllipsePart2();

    TMatrixD *getHomogeneousWDaughterEllipse();
    TMatrixD *getExtendedWDaughterEllipse();

    void getWDaughterMomentum();

    void calcWDaughterEllipse();
    void calcExtendedWDaughterEllipse();

};

