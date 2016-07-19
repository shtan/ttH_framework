#ifndef TOPSYSTEMCHISQUARE
#define TOPSYSTEMCHISQUARE

#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "TString.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/GenVector/LorentzVector.h"

#include "WDaughterEllipseCalculator.h"

#include "commonstruct.h"

using namespace std;
using namespace ROOT::Math;
typedef LorentzVector<PxPyPzE4D<double>> XYZTLorentzVector;

class topSystemChiSquare
{

    // private:
  protected:
    commonstruct::top_system &topsys;
    // public:

    bool rangeFlag_;

    //moved to public
    WDaughterEllipseCalculator WDaughter2Calc_;
    void resetWDaughter2();


    int& debug_verbosity;

  public:
    topSystemChiSquare( commonstruct::top_system&, int& );
    
    virtual ~topSystemChiSquare();

    //WDaughterEllipseCalculator WDaughter2Calc_;

    void preSetupWDaughter2Ellipse();
    void setupWDaughter2EllipsePart2();
    void calcWDaughter2Ellipse();
    void calc_hperp_nperp();
    void setEllipseAngle();
    double getZ2();

    bool hasTopMassRange() { return rangeFlag_; };

    void calcTopMassRange();


    TMatrixD *getHomogeneousWDaughterEllipse();

    void printTopConstituents();
    void printWDaughter2();
    void printChiSquareInfo();
};


#endif
