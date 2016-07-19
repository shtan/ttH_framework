#include "topSystemChiSquare.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH1D.h"
#include "Math/QuantFuncMathCore.h"

using namespace commonstruct;

topSystemChiSquare::topSystemChiSquare( top_system & topsystem, int &debug )
    : topsys(topsystem), WDaughter2Calc_( topsys, debug ), debug_verbosity(debug)
{
    if (debug_verbosity >= 2)
        cout << "starting topSystemChiSquare constructor" << endl;
    //setupWDaughter2Ellipse();
}

topSystemChiSquare::~topSystemChiSquare() {}


void topSystemChiSquare::setEllipseAngle()
{
   resetWDaughter2();
}

void topSystemChiSquare::resetWDaughter2()
{
    // cout << "Calculating the second W daughter momentum at angle: " << theta
    // << endl;
    WDaughter2Calc_.getWDaughterMomentum();
}

void topSystemChiSquare::preSetupWDaughter2Ellipse()
{
    WDaughter2Calc_.preSetupEllipse();
}

void topSystemChiSquare::setupWDaughter2EllipsePart2()
{
    WDaughter2Calc_.setupEllipsePart2();
}


void topSystemChiSquare::calcWDaughter2Ellipse()
{
    // cout << "Calculating topmassrange" <<endl;
    calcTopMassRange();
    // cout << "Calculating the ellipse in homogenous representation" << endl;
    WDaughter2Calc_.calcWDaughterEllipse();
    // cout << "Calculating the ellipse in extended representation" << endl;
    WDaughter2Calc_.calcExtendedWDaughterEllipse();
}

void topSystemChiSquare::calc_hperp_nperp()
{
    WDaughter2Calc_.calcWDaughterEllipse();
    // cout << "Calculating the ellipse in extended representation" << endl;
    WDaughter2Calc_.calcExtendedWDaughterEllipse();

}

//double topSystemChiSquare::getZ2(double mTop, double mW, double mWDaughter2)
double topSystemChiSquare::getZ2()
{
    //return WDaughter2Calc_.getZ2(mTop, mW, mWDaughter2);
    return WDaughter2Calc_.getZ2();
}

TMatrixD *topSystemChiSquare::getHomogeneousWDaughterEllipse()
{
    return WDaughter2Calc_.getHomogeneousWDaughterEllipse();
}

void topSystemChiSquare::calcTopMassRange()
{
    // Find the ranges of squared top mass values in which Z^2 is positive:
    // Z^2=0 is a polynomial equation of degree 2 in mTop^2

    // First check whether the calculation has been done
    //double currentZ2 = getZ2(topsys.input.mTop_central + topsys.input.mTop_width * topsys.vars.delta_mTop,
    //                         topsys.input.mW_central + topsys.input.mW_width * topsys.vars.delta_mW, topsys.input.Wd2_m);
    double currentZ2 = getZ2();

    if (rangeFlag_ && currentZ2 > 0.) {
        // cout << "Top mass range already calculated" << endl;
        return;
    }

    rangeFlag_ = false;

    // cout << "Using pole mass " << topsys.input.mTop_central << " GeV" << endl;

    double mW = topsys.input.mW_central + topsys.input.mW_width * topsys.vars.delta_mW;
    double mW2 = mW * mW;
    double WDaughter2Mass2 = topsys.input.Wd2_m * topsys.input.Wd2_m;

    double bJetE2 = topsys.calc.b_e() * topsys.calc.b_e();
    double bJetP2 = topsys.calc.b_px() * topsys.calc.b_px() + topsys.calc.b_py() * topsys.calc.b_py() + topsys.calc.b_pz() * topsys.calc.b_pz();
    double bJetP = sqrt(bJetP2);

    double WDaughter1P2 = topsys.calc.Wd1_px() * topsys.calc.Wd1_px() +
                          topsys.calc.Wd1_py() * topsys.calc.Wd1_py() +
                          topsys.calc.Wd1_pz() * topsys.calc.Wd1_pz();
    double WDaughter1P = sqrt(WDaughter1P2);
    double WDaughter1P4 = WDaughter1P2 * WDaughter1P2;
    double WDaughter1E2 = topsys.calc.Wd1_e() * topsys.calc.Wd1_e();
    double WDaughter1E3 = WDaughter1E2 * topsys.calc.Wd1_e();
    double WDaughter1E4 = WDaughter1E2 * WDaughter1E2;

    double c = (topsys.calc.Wd1_px() * topsys.calc.b_px() + topsys.calc.Wd1_py() * topsys.calc.b_py() +
                topsys.calc.Wd1_pz() * topsys.calc.b_pz());
    double s2 = 1. - c * c / (bJetP2 * WDaughter1P2);
    c /= (bJetP * WDaughter1P);

    double mTopEdge1 =
        (topsys.calc.b_e() * WDaughter1E3 -
         topsys.calc.b_e() * topsys.calc.Wd1_e() * (WDaughter2Mass2 - mW2 + WDaughter1P2) +
         WDaughter1E2 * (bJetE2 + mW2 - bJetP * (bJetP + c * WDaughter1P)) +
         WDaughter1P * ((-bJetE2 - mW2 + bJetP2) * WDaughter1P +
                        c * bJetP * (WDaughter2Mass2 - mW2 + WDaughter1P2)) -
         sqrt((WDaughter1E4 + pow(WDaughter2Mass2 - mW2, 2) +
               2 * (WDaughter2Mass2 + mW2) * WDaughter1P2 + WDaughter1P4 -
               2 * WDaughter1E2 * (WDaughter2Mass2 + mW2 + WDaughter1P2)) *
              (pow(c * topsys.calc.Wd1_e() * bJetP - topsys.calc.b_e() * WDaughter1P, 2) +
               bJetP2 * (topsys.calc.Wd1_e() - WDaughter1P) *
                   (topsys.calc.Wd1_e() + WDaughter1P) * s2))) /
        ((topsys.calc.Wd1_e() - WDaughter1P) * (topsys.calc.Wd1_e() + WDaughter1P));
    double mTopEdge2 =
        (topsys.calc.b_e() * WDaughter1E3 -
         topsys.calc.b_e() * topsys.calc.Wd1_e() * (WDaughter2Mass2 - mW2 + WDaughter1P2) +
         WDaughter1E2 * (bJetE2 + mW2 - bJetP * (bJetP + c * WDaughter1P)) +
         WDaughter1P * ((-bJetE2 - mW2 + bJetP2) * WDaughter1P +
                        c * bJetP * (WDaughter2Mass2 - mW2 + WDaughter1P2)) +
         sqrt((WDaughter1E4 + pow(WDaughter2Mass2 - mW2, 2) +
               2 * (WDaughter2Mass2 + mW2) * WDaughter1P2 + WDaughter1P4 -
               2 * WDaughter1E2 * (WDaughter2Mass2 + mW2 + WDaughter1P2)) *
              (pow(c * topsys.calc.Wd1_e() * bJetP - topsys.calc.b_e() * WDaughter1P, 2) +
               bJetP2 * (topsys.calc.Wd1_e() - WDaughter1P) *
                   (topsys.calc.Wd1_e() + WDaughter1P) * s2))) /
        ((topsys.calc.Wd1_e() - WDaughter1P) * (topsys.calc.Wd1_e() + WDaughter1P));

    // assert(mTopEdge1==mTopEdge1);
    // assert(mTopEdge2==mTopEdge2);

    double minRoot = min(mTopEdge1, mTopEdge2);
    double maxRoot = max(mTopEdge1, mTopEdge2);

    mTopEdge1 = minRoot;
    mTopEdge2 = maxRoot;

    //  if(mTopEdge1<=-DBL_MAX)
    if (mTopEdge1 <= 0) {
        cout << "infinite lower edge" << endl;
        // cout << "first W daughter energy is " << topsys.calc.Wd1_e() << endl;
        // cout << "first W daughter momentum is " << WDaughter1P << endl;
        // mTopEdge1 = 0;
        // return;
    }
    //  if(mTopEdge2>= DBL_MAX)
    if (mTopEdge2 >= 500) {
        cout << "infinite upper edge" << endl;
        // cout << "first W daughter energy is " << topsys.calc.Wd1_e() << endl;
        // cout << "first W daughter momentum is " << WDaughter1P << endl;
        // mTopEdge2 = 500;
        // return;
    }

    // double tempTheta=theta_;

    cout << "Printing out some starting values" << endl;

    // cout << "Current top mass is " << tempTopMass << endl;
    //cout << "Current ellipse angle is " << getEllipseAngle() << endl;

    double Z2, mTop;
    topsys.vars.has_high_edge = false;

    if (mTopEdge1 > 0 && mTopEdge2 > 0) {
        // both roots are positive -> 3 top mass intervals in which to check
        // sign of Z^2

        mTopEdge1 = sqrt(mTopEdge1);
        mTopEdge2 = sqrt(mTopEdge2);

        bool pos1 = false, pos2 = false, pos3 = false;

        mTop = 0.5 * mTopEdge1;
        topsys.vars.delta_mTop = mTop - topsys.input.mTop_central;
        Z2 = getZ2();
        if (Z2 > 0)
            pos1 = true;

        mTop = 0.5 * (mTopEdge1 + mTopEdge2);
        topsys.vars.delta_mTop = mTop - topsys.input.mTop_central;
        Z2 = getZ2();
        if (Z2 > 0)
            pos2 = true;

        mTop = 2. * mTopEdge2;
        topsys.vars.delta_mTop = mTop - topsys.input.mTop_central;
        Z2 = getZ2();
        if (Z2 > 0)
            pos3 = true;

        if (pos1 && pos2 && pos3) {
            cout << "Z^2 is positive on all 3 intervals" << endl;
        } else if (pos1 && pos2 && !pos3) {
            cout << "Z^2 is positive on 2 consecutive intervals below "
                 << mTopEdge2 << " and negative above" << endl;
        } else if (pos1 && !pos2 && pos3) {
            // cout << "Z^2 is positive below " << mTopEdge1 << " and above " <<
            // mTopEdge2 << endl;
            if (topsys.input.mTop_central < mTopEdge1) {
                rangeFlag_ = true;
                topsys.vars.mTop_edge_low = 0.;
                topsys.vars.mTop_edge_high = mTopEdge1;
                topsys.vars.has_high_edge = true;
            } else if (topsys.input.mTop_central > mTopEdge1 && topsys.input.mTop_central < mTopEdge2) {
                rangeFlag_ = true;
                if (topsys.input.mTop_central - mTopEdge1 < mTopEdge2 - topsys.input.mTop_central) {
                    topsys.vars.mTop_edge_low = 0.;
                    topsys.vars.mTop_edge_high = mTopEdge1;
                    topsys.vars.has_high_edge = true;
                } else {
                    topsys.vars.mTop_edge_low = mTopEdge2;
                    topsys.vars.mTop_edge_high = -1.;
                    topsys.vars.has_high_edge = false;
                }
            } else if (topsys.input.mTop_central > mTopEdge2) {
                rangeFlag_ = true;
                topsys.vars.mTop_edge_low = mTopEdge2;
                topsys.vars.mTop_edge_high = -1.;
                topsys.vars.has_high_edge = false;
            } else {
                cout << "Logic fail: check the roots" << endl;
            }
        } else if (pos1 && !pos2 && !pos3) {
            cout << "Z^2 is positive below " << mTopEdge1
                 << " and negative on 2 consecutive intervals above" << endl;
        } else if (!pos1 && pos2 && pos3) {
            cout << "Z^2 is negative below " << mTopEdge1
                 << " and positive on 2 consecutive intervals above" << endl;
        } else if (!pos1 && pos2 && !pos3) {
            // cout << "Z^2 is positive between " << mTopEdge1 << " and " <<
            // mTopEdge2 << endl;
            rangeFlag_ = true;
            topsys.vars.mTop_edge_low = mTopEdge1;
            topsys.vars.mTop_edge_high = mTopEdge2;
            topsys.vars.has_high_edge = true;
        } else if (!pos1 && !pos2 && pos3) {
            cout << "Z^2 is negative on 2 consecutive intervals below "
                 << mTopEdge2 << " and positive above" << endl;
        } else if (!pos1 && !pos2 && !pos3) {
            cout << "Z^2 is negative on all 3 intervals" << endl;
        } else {
            cout << "Logic fail" << endl;
        }

    } else if (mTopEdge2 > 0) {
        // only one positive root -> 2 intervals to check: between 0 and
        // sqrt(mTopEdge2), above sqrt(mTopEdge2)

        mTopEdge2 = sqrt(mTopEdge2);

        bool pos1 = false, pos2 = false;

        mTop = 0.5 * mTopEdge2;
        topsys.vars.delta_mTop = mTop - topsys.input.mTop_central;
        Z2 = getZ2();
        if (Z2 > 0)
            pos1 = true;

        mTop = 2. * mTopEdge2;
        topsys.vars.delta_mTop = mTop - topsys.input.mTop_central;
        Z2 = getZ2();
        if (Z2 > 0)
            pos2 = true;

        if (pos1 && pos2) {
            cout << "Z^2 is positive on both sides of the root: " << mTopEdge2
                 << endl;
        } else if (pos1) {
            rangeFlag_ = true;
            topsys.vars.mTop_edge_low = 0.;
            topsys.vars.mTop_edge_high = mTopEdge2;
            topsys.vars.has_high_edge = true;
        } else if (pos2) {
            rangeFlag_ = true;
            topsys.vars.mTop_edge_low = mTopEdge2;
            topsys.vars.mTop_edge_high = -1.;
            topsys.vars.has_high_edge = false;
        } else {
            cout << "Z^2 is negative on both sides of the root: " << mTopEdge2
                 << endl;
        }
    } else {
        // both roots are negative -> 1 interval to check: mTop>0
        mTop = 100;
        topsys.vars.delta_mTop = mTop - topsys.input.mTop_central;
        Z2 = getZ2();
        if (Z2 > 0) {
            rangeFlag_ = true;
            topsys.vars.mTop_edge_low = 0.;
            topsys.vars.mTop_edge_high = -1.;
            topsys.vars.has_high_edge = false;
        } else {
            cout << "Both roots are negative, and Z^2 is negative above 0"
                 << endl;
        }
    }

    cout << "before rangeflag check" << endl;
    if (!rangeFlag_) {
        cout << "No top mass range where Z^2 is positive was found" << endl;
        return;
    }
    cout << "aftercheck" << endl;

    if (topsys.vars.has_high_edge)
        cout << "Allowed top mass range is [ " << topsys.vars.mTop_edge_low << " , "
             << topsys.vars.mTop_edge_high << " ]" << endl;
    else
        cout << "Allowed top mass range is [ " << topsys.vars.mTop_edge_low << " , +inf ["
             << endl;

    if (topsys.vars.has_high_edge && topsys.vars.mTop_edge_high < topsys.vars.mTop_edge_low)
        cout << "Inverted interval" << endl;

    // Now calculate the top mass delta range

    if (topsys.vars.mTop_edge_low < mW) {
        topsys.vars.delta_mTop_range_low = (mW - topsys.input.mTop_central) / topsys.input.mTop_width;
        cout << "first " << topsys.vars.delta_mTop_range_low << endl;
    } else {
        topsys.vars.delta_mTop_range_low =
            (topsys.vars.mTop_edge_low - topsys.input.mTop_central) / topsys.input.mTop_width; // can be positive or negative
        cout << "mtopedgelow = " << topsys.vars.mTop_edge_low << endl;
        cout << "mtop = " << topsys.input.mTop_central << endl;
        cout << "sigmamtop = " << topsys.input.mTop_width << endl;
        cout << "mtopedge1 = " << mTopEdge1 << endl;
        cout << "mtopedge2 = " << mTopEdge2 << endl;
        cout << topsys.vars.delta_mTop_range_low << endl;
    }

    if (topsys.vars.has_high_edge) {
        topsys.vars.delta_mTop_range_high =
            (topsys.vars.mTop_edge_high - topsys.input.mTop_central) / topsys.input.mTop_width; // can be positive or negative
    } else {
        topsys.vars.delta_mTop_range_high = -1.;
    }

    cout << "Setting lower edge delta to " << topsys.vars.delta_mTop_range_low << endl;
    if (topsys.vars.has_high_edge)
        cout << "Setting upper edge delta to " << topsys.vars.delta_mTop_range_high << endl;
    //
    cout << "setting low edge to " << topsys.input.mTop_central + topsys.vars.delta_mTop_range_low * topsys.input.mTop_width
         << endl;
    if (topsys.vars.has_high_edge)
        cout << "setting high edge to "
             << topsys.input.mTop_central + topsys.vars.delta_mTop_range_high * topsys.input.mTop_width << endl;

    // reset to initial mass value

    if (currentZ2 <= 0) {
        cout << "Setting the top mass to a value that makes the current Z^2 "
                "positive"
             << endl;

        //SM's question: isn't the whole point of finding the range that Z2 should be positive within this range?
        //So if we set delta_mTop to be the low point of this range, the resulting Z2 should be positive:
        //shouldn't need to keep adding little bits to delta_mTop in order to achieve positive Z2.
        topsys.vars.delta_mTop = topsys.vars.delta_mTop_range_low;
        currentZ2 = getZ2(); //added by SM
        int whilecount = 0;
        while (currentZ2 <= 0 and whilecount <= 10000) {
            cout << "mtop is " << topsys.input.mTop_central << endl;
            topsys.vars.delta_mTop += 1.e-3 * topsys.input.mTop_central;
            cout << "New deltaMTop is " << topsys.vars.delta_mTop << endl;
            cout << "New top mass is " << topsys.input.mTop_central + topsys.input.mTop_width * topsys.vars.delta_mTop
                 << endl;
            //setupWDaughter2EllipsePart2();// already called in getZ2()
            // WDaughter2Calc_.setupEllipse(topsys.input.mTop_central+topsys.input.mTop_width*topsys.vars.delta_mTop,mW,topsys.input.Wd2_m);
            currentZ2 =
                getZ2();
            cout << "Current Z^2 is " << currentZ2 << endl;
            whilecount++;
        }
        // WDaughter2Calc_.calcWDaughterEllipse();
        // WDaughter2Calc_.calcExtendedWDaughterEllipse();
    }

    else {
        cout << "Resetting to starting values" << endl;

        topsys.vars.delta_mTop = 0.; //added by SM.  SM's note: not sure about this?
        setupWDaughter2EllipsePart2();
        // WDaughter2Calc_.setupEllipse(topsys.input.mTop_central+topsys.input.mTop_width*topsys.vars.delta_mTop,mW,topsys.input.Wd2_m);
        // calcWDaughter2Ellipse();
        // WDaughter2Calc_.calcWDaughterEllipse();
        // WDaughter2Calc_.calcExtendedWDaughterEllipse();
        cout << "Resetting second W daughter" << endl;
        // if(tempTheta!=0) resetWDaughter2(tempTheta);

        // printWDaughter2();

        // cout << "Re-printing values -- they should match the starting values
        // printed above" << endl;

        // cout << "Current top mass is " << WDaughter2Calc_.getTopMass() <<
        // endl;

        // if(tempTheta!=getEllipseAngle()) cout << "Issue resetting ellipse
        // angle to starting value" << endl;

        // cout << "Current ellipse angle is " << tempTheta << endl;
    }
}


void topSystemChiSquare::printChiSquareInfo()
{
    cout << "PRINTING CHI2 INFO" << endl;

    cout << "leptonic = " << topsys.input.leptonic << endl;
    cout
        << "b-jet current: "
        << "\npt = " << topsys.calc.b_pt() << "\neta = " << topsys.calc.b_eta() << "\nphi = " << topsys.calc.b_phi()
        << "\nm  = " << topsys.calc.b_m()
        << endl;

    cout << "b-jet inputs: "
        << "\npt = " << topsys.input.b_pt << "\neta = " << topsys.input.b_eta << "\nphi = " << topsys.input.b_phi
        << "\nm  = " << topsys.input.b_m
        << endl;

    cout << "b-jet deltas: "
        << "\npt = " << topsys.vars.b_delta_pt << "\neta = " << topsys.vars.b_delta_eta << "\nphi = " << topsys.vars.b_delta_phi
        << endl;

    cout << "b-jet widths: "
        << "\npt = " << topsys.input.b_pt_width << "\neta = " << topsys.input.b_eta_width << "\nphi = " << topsys.input.b_phi_width
        << endl;

    cout << "b-jet chi2 = " << topsys.calc.b_chi2() << endl;

/*    cout
        << "first W daughter: "
        << "\npx = " << topsys.calc.Wd1_px() << "\npy = " << topsys.calc.Wd1_py()
        << "\npz = " << topsys.calc.Wd1_pz() << "\ne  = " << topsys.calc.Wd1_e()
        << endl;*/

}

void topSystemChiSquare::printTopConstituents()
{
    /*cout << "lepton flag = " << topsys.input.leptonic << endl;
    cout << "Top decay products:" << endl;
    cout << "input b pt = " << topsys.input.b_pt<< endl;
    topsys.vars.b_delta_pt = 23.3;
    cout<< "delta b pt = " << topsys.vars.b_delta_pt << endl;
    //cout << "testfunc " << topsys.input.testfunc() << endl;
    cout << "test " << topsys.testt.b_ptt() << endl;
    cout << "width " << topsys.input.b_pt_width << endl;
    cout <<"calc " << topsys.calc.vars.b_delta_pt << endl;
    cout << "b pt = " << topsys.calc.b_pt()<< endl;*/
    //  cout << "In the derived class:" << endl;
    cout << endl;
    cout << "PRINTING TOP CONSTITUENTS" << endl;
    cout << "leptonic = " << topsys.input.leptonic << endl;
    cout
        << "b-jet: "
        << "\npx = " << topsys.calc.b_px() << "\npy = " << topsys.calc.b_py() << "\npz = " << topsys.calc.b_pz()
        << "\ne  = " << topsys.calc.b_e()
        //       << "\nm  = " <<
        //       sqrt(max(0.,topsys.calc.b_e()*topsys.calc.b_e()-topsys.calc.b_px()*topsys.calc.b_px()-topsys.calc.b_py()*topsys.calc.b_py()-topsys.calc.b_pz()*topsys.calc.b_pz()))
        << endl;

    cout
        << "first W daughter: "
        << "\npx = " << topsys.calc.Wd1_px() << "\npy = " << topsys.calc.Wd1_py()
        << "\npz = " << topsys.calc.Wd1_pz() << "\ne  = " << topsys.calc.Wd1_e()
        //       << "\nm  = " <<
        //       sqrt(max(0.,topsys.calc.Wd1_e()*topsys.calc.Wd1_e()-topsys.calc.Wd1_px()*topsys.calc.Wd1_px()-topsys.calc.Wd1_py()*topsys.calc.Wd1_py()-topsys.calc.Wd1_pz()*topsys.calc.Wd1_pz()))
        << endl;

/*    cout
        << "second W daughter: "
        << "input values: "
        << "\npt = " << topsys.input.Wd2_pt << "\neta = " << topsys.input.Wd2_eta
        << "\nphi = " << topsys.input.Wd2_phi << "\ne  = " << topsys.calc.Wd2_e()
        << "current values: "
        << "\npx = " << topsys.vars.Wd2_px << "\npy = " << topsys.vars.Wd2_py
        << "\npz = " << topsys.vars.Wd2_pz << "\ne  = " << topsys.calc.Wd2_e()
        //       << "\nm  = " <<
        //       sqrt(max(0.,WDaughter2E_*WDaughter2E_-topsys.calc.Wd2_px()*topsys.calc.Wd2_px()-WDaughter2Py_*WDaughter2Py_-WDaughter2Pz_*WDaughter2Pz_))
        //       << "\nreco m = " << sqrt(reconstructed_WDaughter2Mass2_)
        << endl;*/

    printWDaughter2();

    cout << "total Px for this top = " << topsys.calc.total_px() << endl;
    cout << "total Py for this top = " << topsys.calc.total_py() << endl;

    cout << "input values: " << endl;
    cout << "b_px = " << topsys.input.b_px << endl;
    cout << "Wd1_px = " << topsys.input.Wd1_px << endl;
    cout << "Wd2_px = " << topsys.input.Wd2_px << endl;
    cout << "total px = " << topsys.input.total_px << endl;

    cout << "top mass: " << topsys.calc.mTop() << endl;

    cout << "W mass: " << topsys.calc.mW() << endl;

    //  cout << "reconstructed second light quark: "
    //       << "\npx = " << reconstructed_topsys.calc.Wd2_px()
    //       << "\npy = " << reconstructed_WDaughter2Py_
    //       << "\npz = " << reconstructed_WDaughter2Pz_
    //       << "\ne  = " << reconstructed_WDaughter2E_
    //       << endl;
    //
    //  cout << "Second W daughter momentum:"
    //       << "\npx = " << topSystemChiSquare::topsys.calc.Wd2_px()
    //       << "\npy = " << topSystemChiSquare::WDaughter2Py_
    //       << "\npz = " << topSystemChiSquare::WDaughter2Pz_
    //       << endl;
    //  double low, high;
    //  getTopMassRange(low,high);
}

void topSystemChiSquare::printWDaughter2()
{
    cout << "Current second light quark: "
         << "\npt  = " << topsys.calc.Wd2_pt() << "\nphi = " << topsys.calc.Wd2_phi()
         << "\neta = " << topsys.calc.Wd2_eta() << endl;

    cout << "\npx  = " << topsys.vars.Wd2_px << "\npy = " << topsys.vars.Wd2_py
         << "\npz = " << topsys.vars.Wd2_pz << endl;


    cout << "theta is " << topsys.vars.theta << endl;

    //  cout << "reconstructed second light quark: "
    //       << "\npt  = " << reconstructed_WDaughter2Pt_
    //       << "\nphi = " << reconstructed_WDaughter2Phi_
    //       << "\neta = " << reconstructed_WDaughter2Eta_
    //       << endl;
    //
    cout << "dif pt  = " << topsys.calc.Wd2_dif_pt()
         << "\ndif phi = " << topsys.calc.Wd2_dif_phi()
         << "\ndif eta = " << topsys.calc.Wd2_dif_eta() << endl;

    // cout << "second W daughter ellipse in homogeneous coordinates:" << endl;
    // TMatrixD* a=getHomogeneousWDaughterEllipse();
    // a->Print();
}


