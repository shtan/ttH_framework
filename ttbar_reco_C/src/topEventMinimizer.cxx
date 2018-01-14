#include "topEventMinimizer.h"

#include <numeric>

#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPDF.h>

using namespace commonstruct;

topEventMinimizer::topEventMinimizer(big_struct & bigstructure, int& debug)
    : bigstruct(bigstructure), nontops( *(bigstruct.nontops_ptr) ),
    nonTopChiSquare_(lightJetChiSquareMinimumSolver(bigstruct, 0, debug, *(bigstruct.nontops_ptr) )),
    debug_verbosity(debug)
{
    //if (debug_verbosity >= 2)
        cout << "Starting topEventMinimizer constructor" << endl;

    maxConsideredChiSquareRoot_ = 30;

    //nonTopChiSquare_.setupEquations();
    Initialize_minimizers(outerMin_, innerMin_);
    //cout << "initialised minimizers" << endl;

    initialize_best_outer_chiSquares();
    //cout << "initialised chis" << endl;
   
}

void topEventMinimizer::Initialize_minimizers(ROOT::Math::Minimizer *&outer,
                                              ROOT::Math::Minimizer *&inner)
{
    //cout << "inside initialize_minimizers" << endl;
    //ROOT::Math::Minimizer *blah = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Minimize");
    //cout << "before root" << endl;
    outer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Minimize");
    //cout << "root" << endl;
    outer->SetMaxFunctionCalls(10000);
    //cout << "max" << endl;
    outer->SetTolerance(0.02);
    //cout << "tolerance" << endl;
    outer->SetPrintLevel(5);
    //cout << "print" << endl;

    //ROOT::Math::MinimizerOptions inner_opt;
    inner = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    inner->SetMaxFunctionCalls(10000);
    inner->SetTolerance(0.001);
    inner->SetPrintLevel(0);

    //inner_opt.SetMinimizerAlgorithm("Simplex");
    //inner->SetOptions(inner_opt);
}

topEventMinimizer::~topEventMinimizer()
{
    // cout << "destructor" << endl;

    delete outerMin_;
    delete innerMin_;
}


void topEventMinimizer::create_tops()
{
    for (int iTop = 0; iTop < bigstruct.n_tops(); ++iTop) {
        topSysChiSqs_.push_back( topSystemChiSquare( (*bigstruct.tops.at(iTop)), debug_verbosity ) );
    }

}


void topEventMinimizer::reset_best_inner_chiSquares(){
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        (*top)->best_inner_params.Wd2_chi2 = 1.e99;
        (*top)->best_inner_params.mTop_chi2 = 1.e99;
    }

    nontops.best_inner_params.chi2 = 1.e99;
}

void topEventMinimizer::initialize_best_outer_chiSquares(){
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        (*top)->best_outer_params.Wd2_chi2 = 1.e99;
        (*top)->best_outer_params.Wd1_chi2 = 1.e99;
        (*top)->best_outer_params.b_chi2 = 1.e99;
        (*top)->best_outer_params.mW_chi2 = 1.e99;
        (*top)->best_outer_params.mTop_chi2 = 1.e99;
    }

    nontops.best_outer_params.chi2 = 1.e99;

}

void topEventMinimizer::printTopConstituents()
{
    int iTop = 0;
    for (auto it = topSysChiSqs_.begin();
         it != topSysChiSqs_.end(); ++it, ++iTop) {
        cout << "This is top number " << iTop + 1 << endl;
        (*it).printTopConstituents();
        (*it).printWDaughter2();
        (*it).printChiSquareInfo();
    }
}

void topEventMinimizer::checkRecoil_after_fit()
{
    cout << "Check Recoil" << endl;
    cout << "recoil px = " << recoil_px(bigstruct) << endl;
    cout << "recoil py = " << recoil_py(bigstruct) << endl;
    double sum_nontop_px_dif = 0;
    double sum_nontop_py_dif = 0;
    for (int i = 0; i < nontops.input.n_objs; ++i){
        sum_nontop_px_dif += nontops.best_outer_params.jet_dif_px.at(i);
        sum_nontop_py_dif += nontops.best_outer_params.jet_dif_py.at(i);
    }
    cout << "sum nontop px difference = " << sum_nontop_px_dif << endl;
    cout << "sum nontop py difference = " << sum_nontop_py_dif << endl;
}

void topEventMinimizer::printNonTopObjects()
{
    cout << "input objects: " << endl;
    for (int iObj = 0; iObj < bigstruct.nontops_ptr->input.n_objs; iObj++) {
        cout << "Light jet " << iObj + 1
             << "\npt = " << nontops.input.jet_pt.at(iObj)
             << "\neta = " << nontops.input.jet_eta.at(iObj)
             << "\nphi = " << nontops.input.jet_phi.at(iObj) << endl;
        cout << "\npx = " << nontops.calc.jet_px_orig.at(iObj)
             << "\npy = " << nontops.calc.jet_py_orig.at(iObj)
             << "\npz = " << nontops.calc.jet_pz_orig.at(iObj) << endl;
    }

    cout << "best outer objects: " << endl;
    nontops.calc.jet_dif_px_given = nontops.best_outer_params.jet_dif_px;
    nontops.calc.jet_dif_py_given = nontops.best_outer_params.jet_dif_py;
    nontops.calc.jet_dif_pz_given = nontops.best_outer_params.jet_dif_pz;
    double sumPx = 0;
    double sumPy = 0;
    for (int iObj = 0; iObj < bigstruct.nontops_ptr->input.n_objs; iObj++) {
        cout << "Light jet " << iObj + 1
             << "\npx = " << nontops.calc.jet_vec_new().at(iObj).Px()
             << "\npy = " << nontops.calc.jet_vec_new().at(iObj).Py()
             << "\npz = " << nontops.calc.jet_vec_new().at(iObj).Pz() << endl;
        sumPx += nontops.calc.jet_vec_new().at(iObj).Px();
        sumPy += nontops.calc.jet_vec_new().at(iObj).Py();
 
    }

    cout << "sum Px and Py of best outer objects: " << endl;
    cout << "sum Px = " << sumPx;
    cout << "sum Py = " << sumPy;


    cout << "recoil px = " << recoil_px(bigstruct) << endl;
    cout << "recoil py = " << recoil_py(bigstruct) << endl;
}


void topEventMinimizer::calcWDaughterEllipses()
{
    int iTop = 0;
    for (auto thisTopChiSquare = topSysChiSqs_.begin();
         thisTopChiSquare != topSysChiSqs_.end();
         thisTopChiSquare++, iTop++) {
        // cout << "Calculating second W daughter ellipse for top number " <<
        // iTop+1 << endl;
        (*thisTopChiSquare).calcWDaughter2Ellipse();
    }
}


void topEventMinimizer::findStartingValues(int nPoints)
{
    //if (debug_verbosity >= 2)
    //    cout << "find Starting Values" << endl;
    // cout << "Determining starting values for the ellipse angles and top mass
    // deltas" << endl;

    calcWDaughterEllipses();

    double startingChi2 = 1.e9;
    double twoPiOverN = 2. * 3.14159265359 / (double)nPoints;
    int whichTop, step;
    double angle;
    vector<double> angles;
/*    // next 3 lines is my addition
    ellipseAngles_.clear();
    ellipseAngles_.push_back(0.);
    ellipseAngles_.push_back(0.);*/

    //calcWDaughterEllipses();

    // cout << "Beginning loop over all possible angle combinations" << endl;

    for (int i = 0; i < pow(nPoints, bigstruct.n_tops() ); i++) {
        whichTop = 1;
        angles.clear();
        for (int j = 0; j < bigstruct.n_tops(); j++) {
            step = (i / whichTop) % nPoints;
            angle = twoPiOverN * step;
            bigstruct.tops.at(j)->vars.theta = angle;
            // cout << "Setting angle for top " << j+1 << " to " << angle <<
            // endl;
            //(topSysChiSqs_.at(j))->setEllipseAngle();
            (topSysChiSqs_.at(j)).setEllipseAngle();
            whichTop *= nPoints;
            angles.push_back(angle);
        }
        //getDxDyFromEllipses();

        //calcHadronicChiSquare(); // W daughter 2 deltas are calculated in
                                 // getDxDyFromEllipses
        //double thisChi2 =
        //    nonTopChiSquare_.getChiSquare() + getHadronicChiSquare();
        nonTopChiSquare_.calcMin();
        double thisChi2 = nontops.best_innermost_params.chi2 + bigstruct.current_total_hadronic_chi2();
        // cout << "This chi2 is: " << thisChi2 << endl;
        if (thisChi2 < startingChi2) {
            startingChi2 = thisChi2;
            //ellipseAngles_ = angles;
            for (int j = 0; j < bigstruct.n_tops(); j++) {
                bigstruct.tops.at(j)->best_inner_params.theta = angles.at(j);
            }
            
        }
    }

    // cout << "The starting value for the inner chi^2 is " << startingChi2 <<
    // endl;

    //SM's note: Following step is unnecessary: topSysChiSq's ellipse angle is re-set
    //by the minimizer anyway
    // set angles corresponding to minimum chi^2
/*    for (int iTop = 0; iTop < nTops_; iTop++) {
        // cout<<ellipseAngles_.at(iTop)<<endl;
        (topSysChiSqs_.at(iTop))
            ->setEllipseAngle(ellipseAngles_.at(iTop));
    }*/

    //SM's note: I think the following two steps are unnecessary,
    //because the values nonTopChi2_ and hadChi2_ are not used again
    //before the minimization occurs
/*    // cout<<"before setup"<<endl;
    setupNonTopChiSquare();
    // cout<<"before calc hadchi"<<endl;
    calcHadronicChiSquare();*/

    // cout << "Corrected non-top objects:" << endl;
    // nonTopChiSquare_.printResults();
    //cout << "end findstartingvalues" << endl;
}

double topEventMinimizer::innerMinimizationOperator(const double *inputDeltas)
{
     //cout<<"at innermin operator"<<endl;
     //cout << "iteration " << iterrr << endl;
     iterrr++;
    // printTopConstituents();
//    vector<double> ellipseAnglesCurrent;
//    vector<double> topMassDeltasCurrent;

    int i = 0;
    for (int j=0; j < bigstruct.n_tops(); ++j){
        bigstruct.tops.at(j)->vars.theta = inputDeltas[i];
        //cout << "inputDeltas[i] at " << i << " = " << inputDeltas[i] << endl;
        ++i;

        bigstruct.tops.at(j)->vars.delta_mTop = inputDeltas[i];
        //cout << "inputDeltas[i] at " << i << " = " << inputDeltas[i] << endl;

        ++i;
    }


    //Calculate Wd2 momenta using ellipse method, given the theta that is being tried by minuit
/*    for (int j=0; j < bigstruct.n_tops(); ++j){
        //(*it)->setupWDaughter2EllipsePart2();
        if (  == true ){
            return bigstruct.last_total_inner_chi2 + 0.1;
        }
        (*it).calc_hperp_nperp();
        (*it).setEllipseAngle();
    }*/

    //Calculate Wd2 momenta using ellipse method, given the theta that is being tried by minuit
    int jhere = 0;
    for (auto it = topSysChiSqs_.begin(); it != topSysChiSqs_.end(); ++it){
        //(*it)->setupWDaughter2EllipsePart2();
        (*it).calc_hperp_nperp();
        (*it).setEllipseAngle();
        if ( bigstruct.tops.at(jhere)->vars.error_flag == true ){
            return bigstruct.last_total_inner_chi2 + 0.1;
        }
        ++jhere;
    }

    //Calculate minimum non-top chisquare configuration for this theta
    //cout << "beginning clacmin nontop" << endl;
    nonTopChiSquare_.calcMin();
    //cout << "finished clacmin nontop" << endl;

    if (bigstruct.current_total_inner_chi2() < bigstruct.current_best_total_inner_chi2() ){
        int itemp = 0;
        for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
            (*top)->best_inner_params.delta_mTop = (*top)->vars.delta_mTop;
            (*top)->best_inner_params.theta = (*top)->vars.theta;
            (*top)->best_inner_params.Wd2_px = (*top)->vars.Wd2_px;
            (*top)->best_inner_params.Wd2_py = (*top)->vars.Wd2_py;
            (*top)->best_inner_params.Wd2_pz = (*top)->vars.Wd2_pz;

            //print
            /*cout << "best_inner_params being set for top " << itemp << endl;
            cout << "theta = " << (*top)->best_inner_params.theta << endl;
            cout << "Wd2_px = " << (*top)->best_inner_params.Wd2_px << endl;
            cout << "Wd2_py = " << (*top)->best_inner_params.Wd2_py << endl;
            cout << "Wd2_pz = " << (*top)->best_inner_params.Wd2_pz << endl;*/

            (*top)->best_inner_params.Wd2_chi2 = (*top)->calc.Wd2_chi2();
            (*top)->best_inner_params.mTop_chi2 = (*top)->calc.mTop_chi2();
            itemp++;
        }
        nontops.best_inner_params.jet_dif_px = nontops.best_innermost_params.jet_dif_px;
        nontops.best_inner_params.jet_dif_py = nontops.best_innermost_params.jet_dif_py;
        nontops.best_inner_params.jet_dif_pz = nontops.best_innermost_params.jet_dif_pz;

        nontops.best_inner_params.chi2 = nontops.best_innermost_params.chi2;
    }

    bigstruct.last_total_inner_chi2 = bigstruct.current_total_inner_chi2();
    //cout << "end innermin operator" << endl;

    //return innerChi2;
    //cout << bigstruct.current_total_inner_chi2() << endl;
    //cout << innerMin_->Edm() << endl;
    return bigstruct.current_total_inner_chi2();
}

void topEventMinimizer::minimizeNonTopChiSquare()
{
     //cout << "Doing inner minimization" << endl;

    //if (debug_verbosity >= 2)
        //cout << "minimizeNonTopChiSquare" << endl;

    //commonstruct::print_top( *(bigstruct.tops.at(0)) );
    //commonstruct::check_met( bigstruct );

    // cout<<"before infunctor"<<endl;
    // Set up the functor
    ROOT::Math::Functor innerFunc(
        this, &topEventMinimizer::innerMinimizationOperator,
        2 * bigstruct.n_tops()); // 1 angle + 1 top mass delta per top system

    // cout<<"before setfunc in"<<endl;
    // Set up the minimization piece:
    innerMin_->SetFunction(innerFunc);

    // cout<<"before find start val"<<endl;
    // Find starting values for the ellipse angles and top mass deltas
    findStartingValues(50);
    // cout<<"after find start val"<<endl;

    // Set up the parameters

    double step_theta = 0.02 * 3.14159265359;
    double step_mTop = 0.1;
    //const double max = 0.1;
    //double max = 3.14159265359;

    //double step_theta = 0.1 * 3.14159265359;
    //double step_mTop = 0.5;

    int iPar = 0;
    for (int iTop = 0; iTop < bigstruct.n_tops(); iTop++) {
        // cout<<"inside itop loop"<<endl;
        // ellipse angle
        TString parName = "theta_";
        parName += iTop;
        innerMin_->SetVariable(iPar, string(parName), bigstruct.tops.at(iTop)->best_inner_params.theta,
                               step_theta);
        //innerMin_->SetLimitedVariable(iPar, string(parName), bigstruct.tops.at(iTop)->best_inner_params.theta,
        //                       step_theta, -max, max);

        iPar += 1;

        // top mass delta
        bool hasHighEdge;
        double deltaMTopRangeLow, deltaMTopRangeHigh;

        hasHighEdge = bigstruct.tops.at(iTop)->vars.has_high_edge;
        deltaMTopRangeLow = bigstruct.tops.at(iTop)->vars.delta_mTop_range_low;
        deltaMTopRangeHigh = bigstruct.tops.at(iTop)->vars.delta_mTop_range_high;
        double middlerange = (deltaMTopRangeLow + deltaMTopRangeHigh)/2;
        /*cout << "deltaMTopRangeLow = " << deltaMTopRangeLow << endl;
        cout << "deltaMTopRangeHigh = " << deltaMTopRangeHigh << endl;
        cout << "middlerange = " << middlerange << endl;*/
        //deltaMTopRangeHigh = min(bigstruct.tops.at(iTop)->vars.delta_mTop_range_high, maxConsideredChiSquareRoot_);

        //test giving it full range of top mass
//        deltaMTopRangeLow = -30.;
//        deltaMTopRangeHigh = 30.;

        //hasHighEdge = 1;
        //deltaMTopRangeLow = -0.1;
        //deltaMTopRangeHigh = 0.1;


        //(topSysChiSqs_.at(iTop))
        //    ->getTopMassDeltaRange(hasHighEdge, deltaMTopRangeLow,
        //                                 deltaMTopRangeHigh);
        //deltaMTopRangeHigh =
        //    min(deltaMTopRangeHigh, maxConsideredChiSquareRoot_);
        parName = "deltaMTop_";
        parName += iTop;

        innerMin_->SetLimitedVariable(iPar, string(parName),
                                           middlerange, step_mTop,
                                           deltaMTopRangeLow, deltaMTopRangeHigh);


        //Test giving it no limits
//        if (iPar == 3){
/*            innerMin_->SetLowerLimitedVariable(iPar, string(parName),
                                               bigstruct.tops.at(iTop)->best_inner_params.delta_mTop, step_mTop,
                                               0.2);*/
/*        innerMin_->SetLowerLimitedVariable(iPar, string(parName),
                                               0.8, step_mTop,
                                               0.2);
        } else {*/


/*        innerMin_->SetVariable(
            iPar, string(parName), bigstruct.tops.at(iTop)->best_inner_params.delta_mTop, step_mTop);*/
       // }

        
        
        
/*        if (hasHighEdge) {
            // cout << "Current top mass delta is " << topMassDeltas_.at(iTop)
            // << endl;
            // cout << "deltaMTop range is " << deltaMTopRangeLow << " to " <<
            // deltaMTopRangeHigh << endl;
            innerMin_->SetLimitedVariable(
                iPar, string(parName), bigstruct.tops.at(iTop)->best_inner_params.delta_mTop, step_mTop,
                deltaMTopRangeLow, deltaMTopRangeHigh);
            //cout << "HAS HIGH EDGE" << endl;
            //cout << "low = " << deltaMTopRangeLow << endl;
            //cout << "high = " << deltaMTopRangeHigh << endl;
        } else {
            // cout << "Current top mass delta is " << topMassDeltas_.at(iTop)
            // << endl;
            // cout << "deltaMTop lower edge is "<< deltaMTopRangeLow << endl;
            innerMin_->SetLowerLimitedVariable(iPar, string(parName),
                                               bigstruct.tops.at(iTop)->best_inner_params.delta_mTop, step_mTop,
                                               deltaMTopRangeLow);
            //innerMin_->SetLimitedVariable(iPar, string(parName),
            //                                   bigstruct.tops.at(iTop)->best_inner_params.delta_mTop, step_mTop,
            //                                   deltaMTopRangeLow, deltaMTopRangeHigh);
            //cout << "HAS NO HIGH EDGE" << endl;
            //cout << "low = " << deltaMTopRangeLow << endl;
        }*/
        iPar += 1;
    }

    // cout << "Starting the inner minimization" << endl;

    //innerMin_->FixVariable(0);
//    innerMin_->FixVariable(1);
    //innerMin_->FixVariable(2);
//    innerMin_->FixVariable(3);

    //cout << "before innerMin Minimize" << endl;
    iterrr = 0;
    innerMin_->Minimize();
    //cout << "end innerMin operator" << endl;


}

double topEventMinimizer::outerMinimizationOperator(const double *inputDeltas)
{
     //std::cout << "at outermin"<<std::endl;
    // printTopConstituents();
    // reset the inner chi^2 minimum for this outer minimizer step

    //if (debug_verbosity >= 2)
        //cout << "outerMinimizationOperator" << endl;


    reset_best_inner_chiSquares();
    

    int i = 0;
    for (int iTop = 0; iTop < bigstruct.n_tops(); ++iTop){
        bigstruct.tops.at(iTop)->vars.b_delta_pt = inputDeltas[i];
        bigstruct.tops.at(iTop)->vars.b_delta_phi = inputDeltas[i+1];
        bigstruct.tops.at(iTop)->vars.b_delta_eta = inputDeltas[i+2];
        bigstruct.tops.at(iTop)->vars.Wd1_delta_pt = inputDeltas[i+3];
        bigstruct.tops.at(iTop)->vars.Wd1_delta_phi = inputDeltas[i+4];
        bigstruct.tops.at(iTop)->vars.Wd1_delta_eta = inputDeltas[i+5];
        bigstruct.tops.at(iTop)->vars.delta_mW = inputDeltas[i+6];
        /*cout << inputDeltas[i] << endl;
        cout << inputDeltas[i+1] << endl;
        cout << inputDeltas[i+2] << endl;
        cout << inputDeltas[i+3] << endl;
        cout << inputDeltas[i+4] << endl;
        cout << inputDeltas[i+5] << endl;
        cout << inputDeltas[i+6] << endl;*/
        i += 7;

        /*cout << "iTop = " << iTop << endl;
        cout << "b_pt, b_phi, b_eta, wd1_pt, wd1_phi, wd1_eta, delta_mW" << endl;
        cout << bigstruct.tops.at(iTop)->vars.b_delta_pt << endl;
        cout << bigstruct.tops.at(iTop)->vars.b_delta_phi << endl;
        cout << bigstruct.tops.at(iTop)->vars.b_delta_eta << endl;
        cout << bigstruct.tops.at(iTop)->vars.Wd1_delta_pt << endl;
        cout << bigstruct.tops.at(iTop)->vars.Wd1_delta_phi << endl;
        cout << bigstruct.tops.at(iTop)->vars.Wd1_delta_eta << endl;
        cout << bigstruct.tops.at(iTop)->vars.delta_mW << endl;*/
        
        (topSysChiSqs_.at(iTop)).preSetupWDaughter2Ellipse();
        bool hasrange = (topSysChiSqs_.at(iTop)).findTopMassRange();
        //If there is no top mass range for which Z2>0, set the chi2 to a value slightly
        //larger than before, and get out so minuit can try the next step
        if (!hasrange) {
            bigstruct.last_total_outer_chi2 += 100;
            cout << "last total outer chi2 + 100 = " << bigstruct.last_total_outer_chi2 << endl;
            return bigstruct.last_total_outer_chi2;
        }

/*        //Test Z2 values
        cout << "This is top " << iTop << endl;
        for (int mt = -30; mt < 30; ++mt){
            bigstruct.tops.at(iTop)->vars.delta_mTop = 0.1*mt;
            double z2 = (topSysChiSqs_.at(iTop)).getZ2();
            cout << 0.1*mt << " mTop = " << bigstruct.tops.at(iTop)->calc.mTop() << "  Z2 = " << z2 << endl;
        } */


    }

    // Calculate the inner piece
    minimizeNonTopChiSquare();

    if (bigstruct.current_total_outer_chi2() < bigstruct.current_best_total_outer_chi2()){
        
        for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){

            //Set best_outer_params for outer loop variables
            (*top)->best_outer_params.b_delta_pt = (*top)->vars.b_delta_pt;
            (*top)->best_outer_params.b_delta_eta = (*top)->vars.b_delta_eta;
            (*top)->best_outer_params.b_delta_phi = (*top)->vars.b_delta_phi;
            (*top)->best_outer_params.Wd1_delta_pt = (*top)->vars.Wd1_delta_pt;
            (*top)->best_outer_params.Wd1_delta_eta = (*top)->vars.Wd1_delta_eta;
            (*top)->best_outer_params.Wd1_delta_phi = (*top)->vars.Wd1_delta_phi;
            (*top)->best_outer_params.delta_mW = (*top)->vars.delta_mW;

            //Set best_outer_params for inner loop variables
            (*top)->best_outer_params.delta_mTop = (*top)->best_inner_params.delta_mTop;
            (*top)->best_outer_params.theta = (*top)->best_inner_params.theta;
            (*top)->best_outer_params.Wd2_px = (*top)->best_inner_params.Wd2_px;
            (*top)->best_outer_params.Wd2_py = (*top)->best_inner_params.Wd2_py;
            (*top)->best_outer_params.Wd2_pz = (*top)->best_inner_params.Wd2_pz;

            //Set best_outer_params for outer loop chi2's
            (*top)->best_outer_params.b_chi2 = (*top)->calc.b_chi2();
            (*top)->best_outer_params.Wd1_chi2 = (*top)->calc.Wd1_chi2();
            (*top)->best_outer_params.mW_chi2 = (*top)->calc.mW_chi2();

            //Set best_outer_params for inner loop chi2's
            (*top)->best_outer_params.Wd2_chi2 = (*top)->best_inner_params.Wd2_chi2;
            (*top)->best_outer_params.mTop_chi2 = (*top)->best_inner_params.mTop_chi2;
        }

        //Set best_outer_params for non-top objects
        nontops.best_outer_params.jet_dif_px = nontops.best_inner_params.jet_dif_px;
        nontops.best_outer_params.jet_dif_py = nontops.best_inner_params.jet_dif_py;
        nontops.best_outer_params.jet_dif_pz = nontops.best_inner_params.jet_dif_pz;

        //Set best_outer_params for non-top chi2
        nontops.best_outer_params.chi2 = nontops.best_inner_params.chi2;

    }

    bigstruct.last_total_outer_chi2 = bigstruct.current_total_outer_chi2();

    //cout << "current total outer chi2 = " << bigstruct.current_total_outer_chi2() << endl;
    //cout << "theta0 = " << bigstruct.tops.at(0)->best_outer_params.theta << endl;
    //cout << "theta1 = " << bigstruct.tops.at(1)->best_outer_params.theta << endl;
    //cout << "deltaMt0 = " << bigstruct.tops.at(0)->best_outer_params.delta_mTop << endl;
    //cout << "deltaMt1 = " << bigstruct.tops.at(1)->best_outer_params.delta_mTop << endl;
    double innerchi2toprint = nontops.best_outer_params.chi2 + bigstruct.tops.at(0)->best_outer_params.Wd2_chi2 + bigstruct.tops.at(1)->best_outer_params.Wd2_chi2 + bigstruct.tops.at(0)->best_outer_params.mTop_chi2 + bigstruct.tops.at(1)->best_outer_params.mTop_chi2;
    //cout << "inner chi2 = " << innerchi2toprint << endl;

    //cout << "end outermin" << endl;
    return bigstruct.current_total_outer_chi2();
}

void topEventMinimizer::minimizeTotalChiSquare()
{
    //if (debug_verbosity >= 2){
        //cout << "minimizeTotalChiSquare" << endl;
        //printTopConstituents();
        //printNonTopObjects();
    //}

    //     std::cout<<"at min"<<std::endl;
    const int nParameters = 7 * bigstruct.n_tops(); // 3 per b-jet + 3 per W daughter 1 + 1
                                        // per W mass = 7 per top system

    //     std::cout<<"before set functor"<<std::endl;
    // Set up the functor
    ROOT::Math::Functor func(
        this, &topEventMinimizer::outerMinimizationOperator, nParameters);

    //     std::cout<<"before setfunc"<<std::endl;
    // Set up the minimization piece:
    outerMin_->SetFunction(func);

    //     std::cout<<"before set min param"<<std::endl;
    // Setup the minimizer parameters

    double step_bPt = 0.1;
    double step_bPhi = 0.1;
    double step_bEta = 0.1;
    double step_Wd1Pt = 0.1;
    double step_Wd1Phi = 0.1;
    double step_Wd1Eta = 0.1;
    double step_mW = 0.1;

/*    double step_bPt = 0.5;
    double step_bPhi = 0.5;
    double step_bEta = 0.5;
    double step_Wd1Pt = 0.5;
    double step_Wd1Phi = 0.5;
    double step_Wd1Eta = 0.5;
    double step_mW = 0.5;*/

    int iPar = 0, iTop = 0;
    const double max = maxConsideredChiSquareRoot_;
    //const double max = 0.001;
/*    for (auto it = topSysChiSqs_.begin();
         it != topSysChiSqs_.end(); ++it, ++iTop) {
        ostringstream convert; // stream used for the (int) conversion
        convert << iTop;
        const string iTop_str = convert.str();
        const string par1 = "bJetPtDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(iPar, par1, bigstruct.tops.at(iTop)->vars.b_delta_pt, step_bPt,
                                      -max, max);
        ++iPar;
        const string par2 = "bJetPhiDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(iPar, par2, bigstruct.tops.at(iTop)->vars.b_delta_phi,
                                      step_bPhi, -max, max);
        ++iPar;
        const string par3 = "bJetEtaDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(iPar, par3, bigstruct.tops.at(iTop)->vars.b_delta_eta,
                                      step_bEta, -max, max);
        ++iPar;
        const string par4 = "WDaughter1PtDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(
            iPar, par4, bigstruct.tops.at(iTop)->vars.Wd1_delta_pt, step_Wd1Pt, -max, max);
        ++iPar;
        const string par5 = "WDaughter1PhiDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(
            iPar, par5, bigstruct.tops.at(iTop)->vars.Wd1_delta_phi, step_Wd1Phi, -max, max);
        ++iPar;
        const string par6 = "WDaughter1EtaDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(
            iPar, par6, bigstruct.tops.at(iTop)->vars.Wd1_delta_eta, step_Wd1Eta, -max, max);
        ++iPar;
        const string par7 = "deltaMW_" + iTop_str;
        outerMin_->SetLimitedVariable(iPar, par7, bigstruct.tops.at(iTop)->vars.delta_mW, step_mW,
                                      -max, max);
        ++iPar;
    }*/

    for (auto it = topSysChiSqs_.begin();
         it != topSysChiSqs_.end(); ++it, ++iTop) {
        ostringstream convert; // stream used for the (int) conversion
        convert << iTop;
        const string iTop_str = convert.str();
        const string par1 = "bJetPtDelta_" + iTop_str;
        outerMin_->SetVariable(iPar, par1, bigstruct.tops.at(iTop)->vars.b_delta_pt, step_bPt);
        ++iPar;
        const string par2 = "bJetPhiDelta_" + iTop_str;
        outerMin_->SetVariable(iPar, par2, bigstruct.tops.at(iTop)->vars.b_delta_phi,
                                      step_bPhi);
        ++iPar;
        const string par3 = "bJetEtaDelta_" + iTop_str;
        outerMin_->SetVariable(iPar, par3, bigstruct.tops.at(iTop)->vars.b_delta_eta,
                                      step_bEta);
        ++iPar;
        const string par4 = "WDaughter1PtDelta_" + iTop_str;
        outerMin_->SetVariable(
            iPar, par4, bigstruct.tops.at(iTop)->vars.Wd1_delta_pt, step_Wd1Pt);
        ++iPar;
        const string par5 = "WDaughter1PhiDelta_" + iTop_str;
        outerMin_->SetVariable(
            iPar, par5, bigstruct.tops.at(iTop)->vars.Wd1_delta_phi, step_Wd1Phi);
        ++iPar;
        const string par6 = "WDaughter1EtaDelta_" + iTop_str;
        outerMin_->SetVariable(
            iPar, par6, bigstruct.tops.at(iTop)->vars.Wd1_delta_eta, step_Wd1Eta);
        ++iPar;
        const string par7 = "deltaMW_" + iTop_str;
        outerMin_->SetVariable(iPar, par7, bigstruct.tops.at(iTop)->vars.delta_mW, step_mW
                                      );
        ++iPar;
    }

    // std::cout<<"before minimize"<<std::endl;
    // cout << "Starting outer minimization" << endl;

/*    //Fix top masses
    innerMin_->FixVariable(1);
    innerMin_->FixVariable(3);*/

/*    //Fix top 1
    innerMin_->FixVariable(0);
    innerMin_->FixVariable(1);
    outerMin_->FixVariable(0);
    outerMin_->FixVariable(1);
    outerMin_->FixVariable(2);
    outerMin_->FixVariable(3);
    outerMin_->FixVariable(4);
    outerMin_->FixVariable(5);
    outerMin_->FixVariable(6);*/

/*    //Fix top 2
    innerMin_->FixVariable(2);
    innerMin_->FixVariable(3);
    outerMin_->FixVariable(7);
    outerMin_->FixVariable(8);
    outerMin_->FixVariable(9);
    outerMin_->FixVariable(10);
    outerMin_->FixVariable(11);
    outerMin_->FixVariable(12);
    outerMin_->FixVariable(13);*/

/*    //Fix top 1
    //innerMin_->FixVariable(0);
    //innerMin_->FixVariable(1);
    outerMin_->FixVariable(0);
    outerMin_->FixVariable(1);
    outerMin_->FixVariable(2);
    outerMin_->FixVariable(3);
    outerMin_->FixVariable(4);
    outerMin_->FixVariable(5);
    outerMin_->FixVariable(6);

    //Fix top 2
    //innerMin_->FixVariable(2);
    //innerMin_->FixVariable(3);
    outerMin_->FixVariable(7);
    outerMin_->FixVariable(8);
    outerMin_->FixVariable(9);
    outerMin_->FixVariable(10);
    outerMin_->FixVariable(11);
    outerMin_->FixVariable(12);
    outerMin_->FixVariable(13);*/

    //innerMin_->Minimize();
    
    outerMin_->Minimize();

    double* x;
    double* y;
    unsigned int n = 400;
    cout << "before scan" << endl;
    //outerMin_->Scan(6, n, x, y, 60, 100);
/*
    cout << "after scan" << endl;
    TCanvas* canvas = new TCanvas("canvas");
    TGraph *g = new TGraph(400,x,y);

    cout << "after tgraph" << endl;
    g->Draw("ap");
    cout << "after draw" << endl;
    canvas->Print("plotplot.pdf");*/

    //innerMin_->ReleaseVariable(1);
    //innerMin_->ReleaseVariable(3);
    
    //outerMin_->Minimize();

    // std::cout<<"after minimise"<<std::endl;

     cout << "Printing outer min results" << endl;
    outerMin_->SetPrintLevel(1);
    outerMin_->PrintResults();
    cout << "Outer min status is " << outerMin_->Status() << endl;

    // cout << "Printing inner min results" << endl;
    innerMin_->SetPrintLevel(4);
    innerMin_->PrintResults();
    cout << "Inner min status is " << innerMin_->Status() << endl;

    bigstruct.outerMinStatus = outerMin_->Status();
    bigstruct.innerMinStatus = innerMin_->Status();
    bigstruct.outerMin_Edm = outerMin_->Edm();
    bigstruct.innerMin_Edm = innerMin_->Edm();

    setBestValues();

    //plot_chis();

    cout << "THIS IS TOP 1" << endl;
    commonstruct::print_top_all( *(bigstruct.tops.at(0)) );
    cout << "THIS IS TOP 2" << endl;
    commonstruct::print_top_all( *(bigstruct.tops.at(1)) );
    commonstruct::check_met( bigstruct );
 

}

void topEventMinimizer::setBestValues()
{
     //cout << "setBestValues begin" << endl;
    // found: " << chi2Best_ << endl;

    //Sets vars sub-structure elements (i.e. Minuit input parameters)
    //in top_system and nontop_system structures in bigstruct,
    //so that their calculator functions can be used to calculate best-fit values of p, etc.

    cout << "0" << endl;

    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){

        cout << "0.5" << endl;

        (*top)->vars.b_delta_pt = (*top)->best_outer_params.b_delta_pt;

        cout << "0.6" << endl;
        (*top)->vars.b_delta_eta = (*top)->best_outer_params.b_delta_eta;
        (*top)->vars.b_delta_phi = (*top)->best_outer_params.b_delta_phi;
        (*top)->vars.Wd1_delta_pt = (*top)->best_outer_params.Wd1_delta_pt;
        (*top)->vars.Wd1_delta_eta = (*top)->best_outer_params.Wd1_delta_eta;
        (*top)->vars.Wd1_delta_phi = (*top)->best_outer_params.Wd1_delta_phi;

        (*top)->vars.theta = (*top)->best_outer_params.theta;
        (*top)->vars.Wd2_px = (*top)->best_outer_params.Wd2_px;
        (*top)->vars.Wd2_py = (*top)->best_outer_params.Wd2_py;
        (*top)->vars.Wd2_pz = (*top)->best_outer_params.Wd2_pz;

/*        cout << "theta = " << (*top)->vars.theta << endl;
        cout << "Wd2_px = " << (*top)->vars.Wd2_px << endl;
        cout << "Wd2_py = " << (*top)->vars.Wd2_py << endl;
        cout << "Wd2_pz = " << (*top)->vars.Wd2_pz << endl;*/
        cout << "B chi2 from inside = " << (*top)->best_outer_params.b_chi2 << endl;
        cout << "B delta pt = " << (*top)->best_outer_params.b_delta_pt << endl;
        cout << "B delta phi = " << (*top)->best_outer_params.b_delta_phi<< endl;
        cout << "B delta eta = " << (*top)->best_outer_params.b_delta_eta << endl;
        cout << "mW chi2 from inside = " << (*top)->best_outer_params.mW_chi2 << endl;
        cout << "mTop chi2 from inside = " << (*top)->best_outer_params.mTop_chi2 << endl;

        (*top)->vars.delta_mTop = (*top)->best_outer_params.delta_mTop;
        (*top)->vars.delta_mW = (*top)->best_outer_params.delta_mW;
 
    }
 
    cout << "1" << endl;

    nontops.calc.jet_dif_px_given = nontops.best_outer_params.jet_dif_px;
    nontops.calc.jet_dif_py_given = nontops.best_outer_params.jet_dif_py;
    nontops.calc.jet_dif_pz_given = nontops.best_outer_params.jet_dif_pz;

    cout << "2" << endl;

    printTopConstituents();
    cout << "total top Px = " << total_top_px(bigstruct) << endl;
    cout << "total top Py = " << total_top_py(bigstruct) << endl;

    printNonTopObjects();
    checkRecoil_after_fit();

}

TLorentzVector topEventMinimizer::get_b(int iTop)
{
    double pt = bigstruct.tops.at(iTop)->calc.b_pt();
    double eta = bigstruct.tops.at(iTop)->calc.b_eta();
    double phi = bigstruct.tops.at(iTop)->calc.b_phi();
    double m = bigstruct.tops.at(iTop)->calc.b_m();

    return lorentz_maker_pol(pt, eta, phi, m);
}

TLorentzVector topEventMinimizer::get_Wd1(int iTop)
{
    double pt = bigstruct.tops.at(iTop)->calc.Wd1_pt();
    double eta = bigstruct.tops.at(iTop)->calc.Wd1_eta();
    double phi = bigstruct.tops.at(iTop)->calc.Wd1_phi();
    double m = bigstruct.tops.at(iTop)->calc.Wd1_m();

    return lorentz_maker_pol(pt, eta, phi, m);
}

TLorentzVector topEventMinimizer::get_Wd2(int iTop)
{
    double pt = bigstruct.tops.at(iTop)->calc.Wd2_pt();
    double eta = bigstruct.tops.at(iTop)->calc.Wd2_eta();
    double phi = bigstruct.tops.at(iTop)->calc.Wd2_phi();
    double m = bigstruct.tops.at(iTop)->calc.Wd2_m();

    return lorentz_maker_pol(pt, eta, phi, m);
}

TLorentzVector topEventMinimizer::get_W(int iTop)
{
    TLorentzVector W = get_Wd1(iTop) + get_Wd2(iTop);
    return W;
}

TLorentzVector topEventMinimizer::get_top(int iTop)
{
    TLorentzVector t = get_W(iTop) + get_b(iTop);
    return t;
}

TLorentzVector topEventMinimizer::get_nontop_object(int iObj)
{
    return bigstruct.nontops_ptr->calc.jet_vec_new().at(iObj);
}

//Only works if we know the first two non-top objects are the b's from the Higgs.
TLorentzVector topEventMinimizer::get_higgs()
{
    if (bigstruct.nontops_ptr->input.n_objs >= 2){
        TLorentzVector h = get_nontop_object(0) + get_nontop_object(1);
        return h;
    } else {
        cout << "Fewer than 2 non-top objects" << endl;
        TLorentzVector h;
        h.SetPxPyPzE(0,0,0,0);
        return h;
    }
}

double topEventMinimizer::get_best_mW_chi2(int iTop)
{
    double chi2 = bigstruct.tops.at(iTop)->best_outer_params.mW_chi2;
    return chi2;
}

double topEventMinimizer::get_best_mTop_chi2(int iTop)
{
    double chi2 = bigstruct.tops.at(iTop)->best_outer_params.mTop_chi2;
    return chi2;
}

double topEventMinimizer::get_best_b_chi2(int iTop)
{
    double chi2 = bigstruct.tops.at(iTop)->best_outer_params.b_chi2;
    return chi2;
}

double topEventMinimizer::get_best_Wd1_chi2(int iTop)
{
    double chi2 = bigstruct.tops.at(iTop)->best_outer_params.Wd1_chi2;
    return chi2;
}

double topEventMinimizer::get_best_Wd2_chi2(int iTop)
{
    double chi2 = bigstruct.tops.at(iTop)->best_outer_params.Wd2_chi2;
    return chi2;
}

double topEventMinimizer::get_best_nontop_chi2()
{
    double chi2 = bigstruct.nontops_ptr->best_outer_params.chi2;
    return chi2;
}

double topEventMinimizer::get_best_total_had_chi2()
{
    double chi2 = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        chi2 += (*top)->best_outer_params.Wd2_chi2;
    }
    return chi2;
}

double topEventMinimizer::get_best_total_mTop_chi2()
{
    double chi2 = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        chi2 += (*top)->best_outer_params.mTop_chi2;
    }
    return chi2;
}

double topEventMinimizer::get_best_total_topsys_chi2()
{
    double chi2 = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        chi2 += (*top)->best_outer_params.Wd1_chi2;
        chi2 += (*top)->best_outer_params.b_chi2;
        chi2 += (*top)->best_outer_params.mW_chi2;
    }
    return chi2;
}

double topEventMinimizer::get_best_total_chi2()
{
    double chi2 = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        chi2 += (*top)->best_outer_params.Wd1_chi2;
        chi2 += (*top)->best_outer_params.b_chi2;
        chi2 += (*top)->best_outer_params.mW_chi2;
        chi2 += (*top)->best_outer_params.mTop_chi2;
        chi2 += (*top)->best_outer_params.Wd2_chi2;
    }
    chi2 += bigstruct.nontops_ptr->best_outer_params.chi2;

    return chi2;
}

int topEventMinimizer::get_inner_min_status()
{
    return bigstruct.innerMinStatus;
}

int topEventMinimizer::get_outer_min_status()
{
    return bigstruct.outerMinStatus;
}

double topEventMinimizer::get_inner_edm()
{
    return bigstruct.innerMin_Edm;
}
double topEventMinimizer::get_outer_edm()
{
    return bigstruct.outerMin_Edm;
}

void topEventMinimizer::plot_chi2s()
{
    plot_chi2_outer(0, bigstruct.tops.at(0)->vars.b_delta_pt, bigstruct.tops.at(0)->input.b_pt_width,
                    bigstruct.tops.at(0)->input.b_pt, bigstruct.tops.at(0)->best_outer_params.b_delta_pt,
                    "Bottom_1_Pt");
    plot_chi2_outer(0, bigstruct.tops.at(0)->vars.b_delta_phi, bigstruct.tops.at(0)->input.b_phi_width,
                    bigstruct.tops.at(0)->input.b_phi, bigstruct.tops.at(0)->best_outer_params.b_delta_phi,
                    "Bottom_1_Phi");
    plot_chi2_outer(0, bigstruct.tops.at(0)->vars.b_delta_eta, bigstruct.tops.at(0)->input.b_eta_width,
                    bigstruct.tops.at(0)->input.b_eta, bigstruct.tops.at(0)->best_outer_params.b_delta_eta,
                    "Bottom_1_Eta");

    plot_chi2_outer(0, bigstruct.tops.at(0)->vars.Wd1_delta_pt, bigstruct.tops.at(0)->input.Wd1_pt_width,
                    bigstruct.tops.at(0)->input.Wd1_pt, bigstruct.tops.at(0)->best_outer_params.Wd1_delta_pt,
                    "Wd11_Pt");
    plot_chi2_outer(0, bigstruct.tops.at(0)->vars.Wd1_delta_phi, bigstruct.tops.at(0)->input.Wd1_phi_width,
                    bigstruct.tops.at(0)->input.Wd1_phi, bigstruct.tops.at(0)->best_outer_params.Wd1_delta_phi,
                    "Wd11_Phi");
    plot_chi2_outer(0, bigstruct.tops.at(0)->vars.Wd1_delta_eta, bigstruct.tops.at(0)->input.Wd1_eta_width,
                    bigstruct.tops.at(0)->input.Wd1_eta, bigstruct.tops.at(0)->best_outer_params.Wd1_delta_eta,
                    "Wd11_Eta");

    plot_chi2_outer(1, bigstruct.tops.at(1)->vars.b_delta_pt, bigstruct.tops.at(1)->input.b_pt_width,
                    bigstruct.tops.at(1)->input.b_pt, bigstruct.tops.at(1)->best_outer_params.b_delta_pt,
                    "Bottom_2_Pt");
    plot_chi2_outer(1, bigstruct.tops.at(1)->vars.b_delta_phi, bigstruct.tops.at(1)->input.b_phi_width,
                    bigstruct.tops.at(1)->input.b_phi, bigstruct.tops.at(1)->best_outer_params.b_delta_phi,
                    "Bottom_2_Phi");
    plot_chi2_outer(1, bigstruct.tops.at(1)->vars.b_delta_eta, bigstruct.tops.at(1)->input.b_eta_width,
                    bigstruct.tops.at(1)->input.b_eta, bigstruct.tops.at(1)->best_outer_params.b_delta_eta,
                    "Bottom_2_Eta");

    plot_chi2_outer(1, bigstruct.tops.at(1)->vars.Wd1_delta_pt, bigstruct.tops.at(1)->input.Wd1_pt_width,
                    bigstruct.tops.at(1)->input.Wd1_pt, bigstruct.tops.at(1)->best_outer_params.Wd1_delta_pt,
                    "Wd21_Pt");
    plot_chi2_outer(1, bigstruct.tops.at(1)->vars.Wd1_delta_phi, bigstruct.tops.at(1)->input.Wd1_phi_width,
                    bigstruct.tops.at(1)->input.Wd1_phi, bigstruct.tops.at(1)->best_outer_params.Wd1_delta_phi,
                    "Wd21_Phi");
    plot_chi2_outer(1, bigstruct.tops.at(1)->vars.Wd1_delta_eta, bigstruct.tops.at(1)->input.Wd1_eta_width,
                    bigstruct.tops.at(1)->input.Wd1_eta, bigstruct.tops.at(1)->best_outer_params.Wd1_delta_eta,
                    "Wd21_Eta");


}

void topEventMinimizer::plot_chi2_outer(int iTop, double &delta, const double &width, const double &input, double &best, string name)
{
    int n = 600;
    double spacing = 0.1;
    double min = -30.;
    //double max = 30.;
    //double delt = min;
    double chi2 = 1.e99;
    double varplot[n];
    double chi2plot[n];

    for (int i = 0; i<n; i++){
        delta = i*spacing + min;
        chi2 = bigstruct.current_total_outer_chi2();
        varplot[i] = delta*width + input;
        chi2plot[i] = chi2;
        //delt += spacing;
    }

    setBestValues();

    TCanvas* canvas = new TCanvas("canvas");
    TGraph *g = new TGraph(n,varplot,chi2plot);
    g->SetMarkerStyle(2);

    cout << "drawing " << name << endl;
//    g->Draw("ap");
    g->GetXaxis()->SetTitle(name.c_str());
    g->GetYaxis()->SetTitle("total_chi2");
    g->SetTitle(("Total_chi2 vs " + name).c_str());
    g->Draw("ap");
    cout << "after draw" << endl;
    canvas->Print(("./chi2_plots/chi2_" + name + ".pdf").c_str());
  
}

void topEventMinimizer::plot_Wd22_chi2()
{
    int n = 200;
    double spacing = 0.02;
    double angl1[n];
    double thet1[n];
    double pii = 3.1415926535;
       // = bigstruct.tops.at(0).best_outer_params.theta;

    cout << "theta min is " << bigstruct.tops.at(1)->best_outer_params.theta << endl;

    for (int i = 0; i<n; i++){
        bigstruct.tops.at(1)->vars.theta = i*spacing*pii - pii - pii;
        topSysChiSqs_.at(1).preSetupWDaughter2Ellipse();
        topSysChiSqs_.at(1).setupWDaughter2EllipsePart2();
        topSysChiSqs_.at(1).calc_hperp_nperp();
        topSysChiSqs_.at(1).setEllipseAngle();
        thet1[i] = bigstruct.tops.at(1)->calc.Wd2_chi2();
        angl1[i] = i*spacing*pii - pii - pii;
        cout << "thet1 at i = " << i << " = " << thet1[i] << endl;
    }

    setBestValues();

    TCanvas* canvas = new TCanvas("canvas");
    TGraph *g = new TGraph(n,angl1,thet1);
    g->SetMarkerStyle(2);

    cout << "after tgraph" << endl;
    g->Draw("ap");
     cout << "after draw" << endl;
    canvas->Print("plotplot.pdf");
 
}
