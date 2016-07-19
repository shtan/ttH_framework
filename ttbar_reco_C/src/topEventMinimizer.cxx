#include "topEventMinimizer.h"

#include <numeric>

#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>

using namespace commonstruct;

topEventMinimizer::topEventMinimizer(big_struct & bigstructure, int& debug)
    : bigstruct(bigstructure), nontops( *(bigstruct.nontops_ptr) ),
    nonTopChiSquare_(lightJetChiSquareMinimumSolver(bigstruct, 0, debug, *(bigstruct.nontops_ptr) )),
    debug_verbosity(debug)
{
    if (debug_verbosity >= 2)
        cout << "Starting topEventMinimizer constructor" << endl;

    maxConsideredChiSquareRoot_ = 30;

    //nonTopChiSquare_.setupEquations();
    Initialize_minimizers(outerMin_, innerMin_);

    initialize_best_outer_chiSquares();
   
}

void topEventMinimizer::Initialize_minimizers(ROOT::Math::Minimizer *&outer,
                                              ROOT::Math::Minimizer *&inner)
{
    outer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    outer->SetMaxFunctionCalls(1000000);
    outer->SetTolerance(0.001);
    outer->SetPrintLevel(5);

    inner = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    inner->SetMaxFunctionCalls(1000000);
    inner->SetTolerance(0.001);
    inner->SetPrintLevel(0);
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
    if (debug_verbosity >= 2)
        cout << "find Starting Values" << endl;
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
}

double topEventMinimizer::innerMinimizationOperator(const double *inputDeltas)
{
    // cout<<"at innermin operator"<<endl;
    // printTopConstituents();
//    vector<double> ellipseAnglesCurrent;
//    vector<double> topMassDeltasCurrent;

    int i = 0;
    for (int j=0; j < bigstruct.n_tops(); ++j){
        bigstruct.tops.at(j)->vars.theta = inputDeltas[i];
        ++i;

        bigstruct.tops.at(j)->vars.delta_mTop = inputDeltas[i];
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
        if ( bigstruct.tops.at(jhere)->vars.error_flag == true ){
            return bigstruct.last_total_inner_chi2 + 0.1;
        }
        (*it).calc_hperp_nperp();
        (*it).setEllipseAngle();
        ++jhere;
    }

    //Calculate minimum non-top chisquare configuration for this theta
    nonTopChiSquare_.calcMin();

    if (bigstruct.current_total_inner_chi2() < bigstruct.current_best_total_inner_chi2() ){
        for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
            (*top)->best_inner_params.delta_mTop = (*top)->vars.delta_mTop;
            (*top)->best_inner_params.theta = (*top)->vars.theta;
            (*top)->best_inner_params.Wd2_px = (*top)->vars.Wd2_px;
            (*top)->best_inner_params.Wd2_py = (*top)->vars.Wd2_py;
            (*top)->best_inner_params.Wd2_pz = (*top)->vars.Wd2_pz;

            (*top)->best_inner_params.Wd2_chi2 = (*top)->calc.Wd2_chi2();
            (*top)->best_inner_params.mTop_chi2 = (*top)->calc.mTop_chi2();
        }
        nontops.best_inner_params.jet_dif_px = nontops.best_innermost_params.jet_dif_px;
        nontops.best_inner_params.jet_dif_py = nontops.best_innermost_params.jet_dif_py;
        nontops.best_inner_params.jet_dif_pz = nontops.best_innermost_params.jet_dif_pz;

        nontops.best_inner_params.chi2 = nontops.best_innermost_params.chi2;
    }

    bigstruct.last_total_inner_chi2 = bigstruct.current_total_inner_chi2();

    //return innerChi2;
    return bigstruct.current_total_inner_chi2();
}

void topEventMinimizer::minimizeNonTopChiSquare()
{
    // cout << "Doing inner minimization" << endl;

    if (debug_verbosity >= 2)
        cout << "minimizeNonTopChiSquare" << endl;


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

    int iPar = 0;
    for (int iTop = 0; iTop < bigstruct.n_tops(); iTop++) {
        // cout<<"inside itop loop"<<endl;
        // ellipse angle
        TString parName = "theta_";
        parName += iTop;
        innerMin_->SetVariable(iPar, string(parName), bigstruct.tops.at(iTop)->best_inner_params.theta,
                               0.02 * 3.14159265359);
        iPar += 1;

        // top mass delta
        bool hasHighEdge;
        double deltaMTopRangeLow, deltaMTopRangeHigh;

        hasHighEdge = bigstruct.tops.at(iTop)->vars.has_high_edge;
        deltaMTopRangeLow = bigstruct.tops.at(iTop)->vars.delta_mTop_range_low;
        deltaMTopRangeHigh = min(bigstruct.tops.at(iTop)->vars.delta_mTop_range_high, maxConsideredChiSquareRoot_);

        //test giving it full range of top mass
        deltaMTopRangeLow = -30.;
        deltaMTopRangeHigh = 30.;

        //(topSysChiSqs_.at(iTop))
        //    ->getTopMassDeltaRange(hasHighEdge, deltaMTopRangeLow,
        //                                 deltaMTopRangeHigh);
        //deltaMTopRangeHigh =
        //    min(deltaMTopRangeHigh, maxConsideredChiSquareRoot_);
        parName = "deltaMTop_";
        parName += iTop;
        if (hasHighEdge) {
            // cout << "Current top mass delta is " << topMassDeltas_.at(iTop)
            // << endl;
            // cout << "deltaMTop range is " << deltaMTopRangeLow << " to " <<
            // deltaMTopRangeHigh << endl;
            innerMin_->SetLimitedVariable(
                iPar, string(parName), bigstruct.tops.at(iTop)->best_inner_params.delta_mTop, 0.1,
                deltaMTopRangeLow, deltaMTopRangeHigh);
            //cout << "HAS HIGH EDGE" << endl;
            //cout << "low = " << deltaMTopRangeLow << endl;
            //cout << "high = " << deltaMTopRangeHigh << endl;
        } else {
            // cout << "Current top mass delta is " << topMassDeltas_.at(iTop)
            // << endl;
            // cout << "deltaMTop lower edge is "<< deltaMTopRangeLow << endl;
            innerMin_->SetLowerLimitedVariable(iPar, string(parName),
                                               bigstruct.tops.at(iTop)->best_inner_params.delta_mTop, 0.1,
                                               deltaMTopRangeLow);
            //cout << "HAS NO HIGH EDGE" << endl;
            //cout << "low = " << deltaMTopRangeLow << endl;
        }
        iPar += 1;
    }

    // cout << "Starting the inner minimization" << endl;
    innerMin_->Minimize();


}

double topEventMinimizer::outerMinimizationOperator(const double *inputDeltas)
{
    // std::cout << "at outermin"<<std::endl;
    // printTopConstituents();
    // reset the inner chi^2 minimum for this outer minimizer step

    if (debug_verbosity >= 2)
        cout << "outerMinimizationOperator" << endl;


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
        i += 7;
        
        (topSysChiSqs_.at(iTop)).preSetupWDaughter2Ellipse();

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

    return bigstruct.current_total_outer_chi2();
}

void topEventMinimizer::minimizeTotalChiSquare()
{
    if (debug_verbosity >= 2){
        cout << "minimizeTotalChiSquare" << endl;
        //printTopConstituents();
        //printNonTopObjects();
    }

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

    int iPar = 0, iTop = 0;
    const double max = maxConsideredChiSquareRoot_;
    for (auto it = topSysChiSqs_.begin();
         it != topSysChiSqs_.end(); ++it, ++iTop) {
        ostringstream convert; // stream used for the (int) conversion
        convert << iTop;
        const string iTop_str = convert.str();
        const string par1 = "bJetPtDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(iPar, par1, bigstruct.tops.at(iTop)->vars.b_delta_pt, 0.1,
                                      -max, max);
        ++iPar;
        const string par2 = "bJetPhiDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(iPar, par2, bigstruct.tops.at(iTop)->vars.b_delta_phi,
                                      0.1, -max, max);
        ++iPar;
        const string par3 = "bJetEtaDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(iPar, par3, bigstruct.tops.at(iTop)->vars.b_delta_eta,
                                      0.1, -max, max);
        ++iPar;
        const string par4 = "WDaughter1PtDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(
            iPar, par4, bigstruct.tops.at(iTop)->vars.Wd1_delta_pt, 0.1, -max, max);
        ++iPar;
        const string par5 = "WDaughter1PhiDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(
            iPar, par5, bigstruct.tops.at(iTop)->vars.Wd1_delta_phi, 0.1, -max, max);
        ++iPar;
        const string par6 = "WDaughter1EtaDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(
            iPar, par6, bigstruct.tops.at(iTop)->vars.Wd1_delta_eta, 0.1, -max, max);
        ++iPar;
        const string par7 = "deltaMW_" + iTop_str;
        outerMin_->SetLimitedVariable(iPar, par7, bigstruct.tops.at(iTop)->vars.delta_mW, 0.1,
                                      -max, max);
        ++iPar;
    }

    // std::cout<<"before minimize"<<std::endl;
    // cout << "Starting outer minimization" << endl;
    outerMin_->Minimize();

    // std::cout<<"after minimise"<<std::endl;

    // cout << "Printing outer min results" << endl;
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

    setBestValues();
}

void topEventMinimizer::setBestValues()
{
    // cout << "Setting delta values corresponding to the minimum total chi^2 I
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


