#ifndef commonstruct_h
#define commonstruct_h

#include <vector>
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "commonfunc.h"

using namespace std;
using namespace commonfunc;

namespace commonstruct
{

//struct eventInfo {
    

//};

struct top_system {
/*
    //INPUT -----------------------------------------------------------
    //b jet
    const double b_pt_input, b_eta_input, b_phi_input, b_e_input;
    const double b_pt_width, b_eta_width, b_phi_width, b_e_width;

    //W Daughter 1
    const double Wd1_pt_input, Wd1_eta_input, Wd1_phi_input, Wd1_e_input;
    const double Wd1_pt_width, Wd1_eta_width, Wd1_phi_width, Wd1_e_width;

    //W Daughter 2
    const double Wd2_pt_input, Wd2_eta_input, Wd2_phi_input, Wd2_e_input;
    const double Wd2_pt_width, Wd2_eta_width, Wd2_phi_width, Wd2_e_width;

    //Central value and width for Breit Wigner function for Top and W masses
    const double mTop_central, mW_central, mTop_width, mW_width;
*/

    //VARIABLES THAT ARE CHANGED DURING MINIMISATION ----------------------

    public:

        struct inputs {
            //Leptonic flag
            const bool leptonic;

            //b jet
            const double b_pt, b_eta, b_phi, b_m;
            const double b_pt_width, b_eta_width, b_phi_width;

            //W Daughter 1
            const double Wd1_pt, Wd1_eta, Wd1_phi, Wd1_m;
            const double Wd1_pt_width, Wd1_eta_width, Wd1_phi_width;

            //Central value and width for Breit Wigner function for Top and W masses
            const double mTop_central, mTop_width, mW_central, mW_width;

            //W Daughter 2 (only used for hadronic tops)
            const double Wd2_pt, Wd2_eta, Wd2_phi, Wd2_m;
            const double Wd2_pt_width, Wd2_eta_width, Wd2_phi_width;

            //Convert to Lorentz Vector
            const TLorentzVector b_vec = lorentz_maker_pol(b_pt, b_eta, b_phi, b_m);
            const TLorentzVector Wd1_vec = lorentz_maker_pol(Wd1_pt, Wd1_eta, Wd1_phi, Wd1_m);
            const TLorentzVector Wd2_vec = lorentz_maker_pol(Wd2_pt, Wd2_eta, Wd2_phi, Wd2_m);

            //Convert to cartesian coords
            const double b_px = b_vec.Px();
            const double b_py = b_vec.Py();
            const double b_pz = b_vec.Pz();
            const double b_e = b_vec.E();
            const double Wd1_px = Wd1_vec.Px();
            const double Wd1_py = Wd1_vec.Py();
            const double Wd1_pz = Wd1_vec.Pz();
            const double Wd1_e = Wd1_vec.E();
            const double Wd2_px = Wd2_vec.Px();
            const double Wd2_py = Wd2_vec.Py();
            const double Wd2_pz = Wd2_vec.Pz();
            const double Wd2_e = Wd2_vec.E();

            //Calculate total px, py, pz of the top system for this top
            const double total_px = b_vec.Px() + Wd1_vec.Px() + Wd2_vec.Px();
            const double total_py = b_vec.Py() + Wd1_vec.Py() + Wd2_vec.Py();
            const double total_pz = b_vec.Pz() + Wd1_vec.Pz() + Wd2_vec.Pz();
         

            inputs( const bool &lepbool, 
                    const double &bPt, const double &bEta, const double &bPhi, const double &bM,
                    const double &bPtWidth, const double &bEtaWidth, const double &bPhiWidth,
                    const double &Wd1Pt, const double &Wd1Eta, const double &Wd1Phi, const double &Wd1M,
                    const double &Wd1PtWidth, const double &Wd1EtaWidth, const double &Wd1PhiWidth,
                    const double &mTopCentral, const double &mTopWidth, const double &mWCentral, const double &mWWidth,                                const double &Wd2Pt, const double &Wd2Eta, const double &Wd2Phi, const double &Wd2M,
                    const double &Wd2PtWidth, const double &Wd2EtaWidth, const double &Wd2PhiWidth)
                :   leptonic(lepbool), 
                    b_pt(bPt), b_eta(bEta), b_phi(bPhi), b_m(bM),
                    b_pt_width(bPtWidth), b_eta_width(bEtaWidth), b_phi_width(bPhiWidth),
                    Wd1_pt(Wd1Pt), Wd1_eta(Wd1Eta), Wd1_phi(Wd1Phi), Wd1_m(Wd1M),
                    Wd1_pt_width(Wd1PtWidth), Wd1_eta_width(Wd1EtaWidth), Wd1_phi_width(Wd1PhiWidth),      
                    mTop_central(mTopCentral), mTop_width(mTopWidth), mW_central(mWCentral), mW_width(mWWidth),
                    Wd2_pt(Wd2Pt), Wd2_eta(Wd2Eta), Wd2_phi(Wd2Phi), Wd2_m(Wd2M),
                    Wd2_pt_width(Wd2PtWidth), Wd2_eta_width(Wd2EtaWidth), Wd2_phi_width(Wd2PhiWidth) {}

        };

        inputs input;


        struct variables {
            //inputs &input;

            public:
                //variables(inputs & inpt) : input(inpt) {}
                variables(){}

                //NOTE: deltas for b and Wd1 are given in factors of sigma, i.e. are dimensionless.
                double b_delta_pt = 0;
                double b_delta_eta = 0;
                double b_delta_phi = 0;
                double Wd1_delta_pt = 0;
                double Wd1_delta_eta = 0;
                double Wd1_delta_phi = 0;
                //double Wd2_delta_pt = 0;
                //double Wd2_delta_eta = 0;
                //double Wd2_delta_phi = 0;

                double delta_mTop = 0;
                double delta_mW = 0;

                //Wd2 Ellipse angle
                double theta = 0;
                //Wd2 Momenta, calculated by WDaughterEllipseCalculator given a value of theta
                double Wd2_px = 0;
                double Wd2_py = 0;
                double Wd2_pz = 0;

                //Top mass range variables
                double mTop_edge_high = -1;
                double mTop_edge_low = -1;
                double delta_mTop_range_low = 0;
                double delta_mTop_range_high = 0;
                bool has_high_edge = 0;

                //error flag to indicate if Z^2 < 0
                bool error_flag = false;

        }; //END variables sub-structure

        variables vars;


        struct calculator{
            inputs &input;
            variables &vars;

            public:
                calculator(inputs & inpt, variables & var) : input(inpt), vars(var) {}

                //Functions to calculate current value of momenta
                double b_pt(){
                    //cout << input.b_pt << endl;
                    //cout << vars.b_delta_pt << endl;
                    //cout << input.b_pt_width << endl;
                    return input.b_pt + vars.b_delta_pt*input.b_pt_width; }
                double b_eta(){ return input.b_eta + vars.b_delta_eta*input.b_eta_width; }
                double b_phi(){ return input.b_phi + vars.b_delta_phi*input.b_phi_width; }
                double Wd1_pt(){ return input.Wd1_pt + vars.Wd1_delta_pt*input.Wd1_pt_width; }
                double Wd1_eta(){ return input.Wd1_eta + vars.Wd1_delta_eta*input.Wd1_eta_width; }
                double Wd1_phi(){ return input.Wd1_phi + vars.Wd1_delta_phi*input.Wd1_phi_width; }

                //Functions to calculate current value of top and W masses
                double mTop(){ return input.mTop_central + vars.delta_mTop*input.mTop_width; }
                double mW(){ return input.mW_central + vars.delta_mW*input.mW_width; }


                //Calculate Cartesian momenta of b and Wd1
                double b_px(){ return (lorentz_maker_pol(b_pt(), b_eta(), b_phi(), input.b_m) ).Px(); }
                double b_py(){ return (lorentz_maker_pol(b_pt(), b_eta(), b_phi(), input.b_m) ).Py(); }
                double b_pz(){ return (lorentz_maker_pol(b_pt(), b_eta(), b_phi(), input.b_m) ).Pz(); }
                double b_e(){ return (lorentz_maker_pol(b_pt(), b_eta(), b_phi(), input.b_m) ).E(); }

                double Wd1_px(){ return (lorentz_maker_pol(Wd1_pt(), Wd1_eta(), Wd1_phi(), input.Wd1_m) ).Px(); }
                double Wd1_py(){ return (lorentz_maker_pol(Wd1_pt(), Wd1_eta(), Wd1_phi(), input.Wd1_m) ).Py(); }
                double Wd1_pz(){ return (lorentz_maker_pol(Wd1_pt(), Wd1_eta(), Wd1_phi(), input.Wd1_m) ).Pz(); }
                double Wd1_e(){ return (lorentz_maker_pol(Wd1_pt(), Wd1_eta(), Wd1_phi(), input.Wd1_m) ).E(); }

                //Mass variables
                double b_m(){ return input.b_m; }
                double Wd1_m(){ return input.Wd1_m; }
                double Wd2_m(){ return input.Wd2_m; }
                double b_m2(){ return pow( b_m(), 2 ); }
                double Wd1_m2(){ return pow( Wd1_m(), 2 ); }
                double Wd2_m2(){ return pow( Wd2_m(), 2 ); }
                double mTop2(){ return pow( mTop(), 2 ); }
                double mW2(){ return pow( mW(), 2 ); }


                //Calculate polar momentum of Wd2
                double Wd2_e(){ return sqrt(max(0., pow(input.Wd2_m,2) + pow(vars.Wd2_px, 2) + pow(vars.Wd2_py, 2) + pow(vars.Wd2_pz, 2) )); }
                double Wd2_pt(){ return (lorentz_maker_car(vars.Wd2_px, vars.Wd2_py, vars.Wd2_pz, Wd2_e() ) ).Pt(); }
                double Wd2_eta(){ return (lorentz_maker_car(vars.Wd2_px, vars.Wd2_py, vars.Wd2_pz, Wd2_e() ) ).Eta(); }
                double Wd2_phi(){ return (lorentz_maker_car(vars.Wd2_px, vars.Wd2_py, vars.Wd2_pz, Wd2_e() ) ).Phi(); }
                //double Wd2_e(){ return (lorentz_maker_car(Wd2_px(), Wd2_py(), Wd2_pz(), input.Wd2_m) ).E(); }


                //For Wd2: function to calculate difference between input and current Wd2 value.
                //Only makes sense for hadronic tops.
                //I call this 'dif' instead of 'delta' because for bJet and Wd1, delta is defined as dimensionless.
                //This 'dif' for Wd2 has dimensions of GeV.
                double Wd2_dif_pt(){ return Wd2_pt() - input.Wd2_pt; }
                double Wd2_dif_eta(){ return Wd2_eta() - input.Wd2_eta; }
                double Wd2_dif_phi(){ 
                    double phidif = Wd2_phi() - input.Wd2_phi;
                    while (phidif > 3.14159265359) {
                        phidif -= 2. * 3.14159265359;
                    }
                    while (phidif < -3.14159265359) {
                        phidif += 2. * 3.14159265359;
                    }
                    return phidif;
                }

                //Function to calculate total Px and Py
                double total_px(){ return b_px() + Wd1_px() + vars.Wd2_px; }
                double total_py(){ return b_py() + Wd1_py() + vars.Wd2_py; }
                double total_pz(){ return b_pz() + Wd1_pz() + vars.Wd2_pz; }


                //Calculate chi-squareds
                double b_chi2(){ return vars.b_delta_pt*vars.b_delta_pt
                                        + vars.b_delta_eta*vars.b_delta_eta
                                        + vars.b_delta_phi*vars.b_delta_phi; }
                double Wd1_chi2(){ return vars.Wd1_delta_pt*vars.Wd1_delta_pt
                                        + vars.Wd1_delta_eta*vars.Wd1_delta_eta
                                        + vars.Wd1_delta_phi*vars.Wd1_delta_phi; }

                double Wd2_chi2(){
                    double toreturn;
                    if (input.leptonic == 0){
                        toreturn = pow(Wd2_dif_pt(), 2) / pow(input.Wd2_pt_width, 2)
                                    + pow(Wd2_dif_eta(), 2) / pow(input.Wd2_eta_width, 2)
                                    + pow(Wd2_dif_phi(), 2) / pow(input.Wd2_phi_width, 2);
                    } else {
                        toreturn = 0;
                    }
                    return toreturn;
                }


                double mW_chi2(){ return breitWignerErr(input.mW_central, input.mW_width, vars.delta_mW); }
                double mTop_chi2(){ return breitWignerErr(input.mTop_central, input.mTop_width, vars.delta_mTop); }

                //The following is wrong, because the current values of mTop and Wd2 may simply be 
                //the ones that have last been tried, rather than the best one in the inner loop
                //double topsys_chi2(){ return b_chi2() + Wd1_chi2() + Wd2_chi2() + mW_chi2() + mTop_chi2(); }

        }; //END calculator sub-structure

        calculator calc;

        struct best_inner_parameters {
            best_inner_parameters() {}

            double delta_mTop = 0;

            //Wd2 Ellipse angle
            double theta = 0;
            //Wd2 Momenta, calculated by WDaughterEllipseCalculator given a value of theta
            double Wd2_px = 0;
            double Wd2_py = 0;
            double Wd2_pz = 0;
 
            double Wd2_chi2 = 1.e99;
            double mTop_chi2 = 1.e99;

        }; //END best_inner_parameters definition

        best_inner_parameters best_inner_params;

        //This is not to store the best parameters, rather the parameters from the last loop.
        //Used so that if Z^2 < 0, we can have minuit "skip" this particular set of tried values
        //best_inner_parameters last_inner_params;

        struct best_outer_parameters {
            best_outer_parameters() {}

            double b_delta_pt = 0;
            double b_delta_eta = 0;
            double b_delta_phi = 0;
            double Wd1_delta_pt = 0;
            double Wd1_delta_eta = 0;
            double Wd1_delta_phi = 0;

            double delta_mTop = 0;
            double delta_mW = 0;

            //Wd2 Ellipse angle
            double theta = 0;
            //Wd2 Momenta, calculated by WDaughterEllipseCalculator given a value of theta
            double Wd2_px = 0;
            double Wd2_py = 0;
            double Wd2_pz = 0;
 
            double b_chi2 = 1.e99;
            double Wd1_chi2 = 1.e99;
            double Wd2_chi2 = 1.e99;
            double mW_chi2 = 1.e99;
            double mTop_chi2 = 1.e99;
            //double topsys_chi2 = 1.e99;

        }; //END best_outer_parameters definition

        best_outer_parameters best_outer_params;

        double testvar;

        struct test {
            inputs *inputt;
            double *testvarr;

            public:
                test(inputs & inpt, double & tv) : inputt(&inpt), testvarr(&tv) {}
                double b_delta_pt = 6;
                double b_delta_eta = 0;
                double b_delta_phi = 0;
                double Wd1_delta_pt = 0;
                double Wd1_delta_eta = 0;
                double Wd1_delta_phi = 0;
                double Wd2_delta_pt = 0;
                double Wd2_delta_eta = 0;
                double Wd2_delta_phi = 0;

                double delta_mTop = 0;
                double delta_mW = 0;

                //Ellipse angle
                double theta = 0;

                double b_ptt(){ return b_delta_pt + inputt->b_pt; }
                double b_pttt() { return b_delta_pt + *testvarr; }
        };
        test testt;


        top_system( const bool &lepbool, 
                    const double &bPt, const double &bEta, const double &bPhi, const double &bM,
                    const double &bPtWidth, const double &bEtaWidth, const double &bPhiWidth,
                    const double &Wd1Pt, const double &Wd1Eta, const double &Wd1Phi, const double &Wd1M,
                    const double &Wd1PtWidth, const double &Wd1EtaWidth, const double &Wd1PhiWidth,
                    const double &mTopCentral, const double &mTopWidth, const double &mWCentral, const double &mWWidth,                                const double &Wd2Pt = 0, const double &Wd2Eta = 0, const double &Wd2Phi = 0, const double &Wd2M = 0,
                    const double &Wd2PtWidth = -99, const double &Wd2EtaWidth = -99, const double &Wd2PhiWidth = -99)
                :   input( lepbool, bPt, bEta, bPhi, bM,
                    bPtWidth, bEtaWidth, bPhiWidth,
                    Wd1Pt, Wd1Eta, Wd1Phi, Wd1M,
                    Wd1PtWidth, Wd1EtaWidth, Wd1PhiWidth,
                    mTopCentral, mTopWidth, mWCentral, mWWidth,
                    Wd2Pt, Wd2Eta, Wd2Phi, Wd2M,
                    Wd2PtWidth, Wd2EtaWidth, Wd2PhiWidth), vars(), calc(input, vars), testt(input, testvar),
                    best_inner_params(), best_outer_params() {}

}; //END top_system definition

struct nontop_system{
    public:
        struct inputs{
            const vector< double > jet_pt;
            const vector< double > jet_eta;
            const vector< double > jet_phi;
            const vector< double > jet_m;
            const vector< double > jet_pt_width;
            const vector< double > jet_eta_width;
            const vector< double > jet_phi_width;

            const int n_objs = (int)(jet_pt.size());


            inputs( const vector< double >& jetPt, const vector< double >& jetEta, const vector< double >& jetPhi,
                    const vector< double >& jetM,
                    const vector< double >& jetPtWidth, const vector< double >& jetEtaWidth,
                    const vector< double >& jetPhiWidth)
                :   jet_pt(jetPt), jet_eta(jetEta), jet_phi(jetPhi), jet_m(jetM),
                    jet_pt_width(jetPtWidth), jet_eta_width(jetEtaWidth), jet_phi_width(jetPhiWidth) {
                    
                    }

        }; //END inputs sub-structure definition
        
        inputs input;


        struct best_params{
            inputs &input;
            public:

                best_params(inputs &inpt) : input(inpt) {}

                //Units of GeV
                vector< double > jet_dif_px = vector< double >(input.n_objs, 0.);
                vector< double > jet_dif_py = vector< double >(input.n_objs, 0.);
                vector< double > jet_dif_pz = vector< double >(input.n_objs, 0.);

                double chi2 = 1.e99;
 
        }; //END best_params definition

        //Best parameters from innermost minimisation (using non-top covariance matrix)
        best_params best_innermost_params;
        //Best parameters from inner loop (where theta and mTop are varied)
        best_params best_inner_params;
        //Best parameters from outer loop
        best_params best_outer_params;

        struct calculator{
            inputs &input;

            public:

                //Set these values in order to calculate new jet_P's for these values of the differences
                vector< double > jet_dif_px_given = vector< double >(input.n_objs, 0.);
                vector< double > jet_dif_py_given = vector< double >(input.n_objs, 0.);
                vector< double > jet_dif_pz_given = vector< double >(input.n_objs, 0.);

                //Original (input) values of jet as a TLorentzVector
                vector< TLorentzVector > jet_vec_orig;// = vector< TLorentzVector >(input.n_objs, 0);


                //Original (input) values of px, py, pz, e
                vector< double > jet_px_orig = vector< double >(input.n_objs, 0.);
                vector< double > jet_py_orig = vector< double >(input.n_objs, 0.);
                vector< double > jet_pz_orig = vector< double >(input.n_objs, 0.);
                vector< double > jet_e_orig = vector< double >(input.n_objs, 0.);

                vector< TLorentzVector > jet_vec_new(){
                    vector< double > jet_px_new = vector< double >(input.n_objs, 0.);
                    vector< double > jet_py_new = vector< double >(input.n_objs, 0.);
                    vector< double > jet_pz_new = vector< double >(input.n_objs, 0.);
                    vector< TLorentzVector > vecnew;// = vector< TLorentzVector >(input.n_objs, 0);
                    vecnew.clear();
                    for (int i = 0; i < input.n_objs; ++i){
                        jet_px_new.at(i) = jet_px_orig.at(i) + jet_dif_px_given.at(i);
                        jet_py_new.at(i) = jet_py_orig.at(i) + jet_dif_py_given.at(i);
                        jet_pz_new.at(i) = jet_pz_orig.at(i) + jet_dif_pz_given.at(i);

                        vecnew.push_back( lorentz_maker_carM( jet_px_new.at(i), jet_py_new.at(i),
                                                            jet_pz_new.at(i), input.jet_m.at(i) ) );

                    }
                    return vecnew;
                }


                calculator(inputs &inpt) : input(inpt) {
                    jet_vec_orig.clear();
                    for (int i = 0; i < input.n_objs; ++i){
                    jet_vec_orig.push_back( lorentz_maker_pol( input.jet_pt.at(i), input.jet_eta.at(i),
                                                            input.jet_phi.at(i), input.jet_m.at(i)) );
                    }
                    for (int i = 0; i < input.n_objs; ++i){
                        jet_px_orig.at(i) = jet_vec_orig.at(i).Px();
                        jet_py_orig.at(i) = jet_vec_orig.at(i).Py();
                        jet_pz_orig.at(i) = jet_vec_orig.at(i).Pz();
                        jet_e_orig.at(i) = jet_vec_orig.at(i).E();
                    }
 
                }

        }; //END calculator sub-structure definition

        calculator calc;

        nontop_system( const vector< double >& jetPt, const vector< double >& jetEta, const vector< double >& jetPhi,
                       const vector< double >& jetM,
                       const vector< double >& jetPtWidth, const vector< double >& jetEtaWidth,
                       const vector< double >& jetPhiWidth )
                   :   input(jetPt, jetEta, jetPhi, jetM, jetPtWidth, jetEtaWidth, jetPhiWidth),
                           best_innermost_params(input), best_inner_params(input),
                           best_outer_params(input), calc(input) {}
 

}; //END nontop_system definition

struct big_struct{
    public:

        big_struct() {}

        vector< top_system* > tops;

        int n_tops(){ return (int)(tops.size()); }

        //Declare a nontop_system object, putting in dummy parameters for now.
        //Will be replaced with an object with the real parameters 
        //when non-top objects are added in topReconstructionFromLHE_core.
        //const vector <double> vec;

        nontop_system* nontops_ptr = NULL;

        //Only one nontop_system object is needed, so this vector should have size 1.
        //The vector is created just as a trick to allow bigstruct to include a nontop_system object,
        //without initializing the nontop_system object right when bigstruct is initialized.
        vector< nontop_system > nontopsvec;

        double MET_px = 0;
        double MET_py = 0;

        double last_total_inner_chi2;
        double last_total_outer_chi2 = 1e20;

        double current_total_hadronic_chi2() {
            double hadchi2 = 0;
            for (auto top = tops.begin(); top != tops.end(); ++top){
                hadchi2 += (*top)->calc.Wd2_chi2();
            }
            return hadchi2;
        }

        double current_total_inner_chi2() {
            double topinnerchi2 = 0;
            for (auto top = tops.begin(); top != tops.end(); ++top){
                topinnerchi2 += (*top)->calc.mTop_chi2();
                topinnerchi2 += (*top)->calc.Wd2_chi2();
            }
            return topinnerchi2 + nontops_ptr->best_innermost_params.chi2;
        }

        double current_best_total_inner_chi2() {
            double besttopinnerchi2 = 0;
            for (auto top = tops.begin(); top != tops.end(); ++top){
                besttopinnerchi2 += (*top)->best_inner_params.mTop_chi2;
                besttopinnerchi2 += (*top)->best_inner_params.Wd2_chi2;
            }
            return besttopinnerchi2 + nontops_ptr->best_inner_params.chi2;
        }

        double current_total_outer_chi2() {
            double topouterchi2 = 0;
            for (auto top = tops.begin(); top != tops.end(); ++top){
                topouterchi2 += (*top)->calc.mW_chi2();
                topouterchi2 += (*top)->calc.Wd1_chi2();
                topouterchi2 += (*top)->calc.b_chi2();
            }
            return topouterchi2 + current_best_total_inner_chi2();
        }
        
        double current_best_total_outer_chi2() {
            double besttopouterchi2 = 0;
            for (auto top = tops.begin(); top != tops.end(); ++top){
                besttopouterchi2 += (*top)->best_outer_params.mTop_chi2;
                besttopouterchi2 += (*top)->best_outer_params.Wd2_chi2;
                besttopouterchi2 += (*top)->best_outer_params.Wd1_chi2;
                besttopouterchi2 += (*top)->best_outer_params.mW_chi2;
                besttopouterchi2 += (*top)->best_outer_params.b_chi2;
            }
            return besttopouterchi2 + nontops_ptr->best_outer_params.chi2;
        }


        int innerMinStatus = -1;
        int outerMinStatus = -1;
        double outerMin_Edm = -1;
        double innerMin_Edm = -1;
        //big_struct() {}

}; //END big_struct definition


inline double adder(top_system & cha){
    TLorentzVector bob;
    bob.SetPtEtaPhiM(cha.input.b_pt, cha.input.b_eta, cha.input.b_phi, cha.input.b_m);
    return bob.Px();
}

inline void check_met(big_struct & bigstruct){
    double totalpx_orig = 0;
    double totalpy_orig = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        totalpx_orig += (*top)->input.total_px;
        totalpy_orig += (*top)->input.total_py;
    }
    cout << "total top input px = " << totalpx_orig << endl;
    cout << "total top input py = " << totalpy_orig << endl;

    double sum_nontop_px = 0;
    double sum_nontop_py = 0;
    for (int i = 0; i < bigstruct.nontops_ptr->input.n_objs; ++i){
        sum_nontop_px += bigstruct.nontops_ptr->calc.jet_px_orig.at(i);
        sum_nontop_py += bigstruct.nontops_ptr->calc.jet_py_orig.at(i);
    }
    cout << "total nontop input px = " << sum_nontop_px << endl;
    cout << "total nontop input py = " << sum_nontop_py << endl;


    cout << "input met px = " << bigstruct.MET_px << endl;
    cout << "input met py = " << bigstruct.MET_py << endl;

    return;
}


inline double recoil_px(big_struct & bigstruct){
    double totalpx = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        totalpx += (*top)->calc.total_px();
    }
    double totalpx_orig = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        totalpx_orig += (*top)->input.total_px;
    }
    totalpx_orig += bigstruct.MET_px;

    return -(totalpx - totalpx_orig);
}

inline double recoil_py(big_struct & bigstruct){
    double totalpy = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        totalpy += (*top)->calc.total_py();
    }
    double totalpy_orig = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        totalpy_orig += (*top)->input.total_py;
    }
    totalpy_orig += bigstruct.MET_py;

    return -(totalpy - totalpy_orig);
}

inline double recoil_pz(big_struct & bigstruct){
    double totalpz = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        totalpz += (*top)->calc.total_pz();
    }
    double totalpz_orig = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        totalpz_orig += (*top)->input.total_pz;
    }

    return -(totalpz - totalpz_orig);
}

inline double total_top_px(big_struct & bigstruct){
    double totalpx = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        totalpx += (*top)->calc.total_px();
    }
    return totalpx;
}

inline double total_top_py(big_struct & bigstruct){
    double totalpy = 0;
    for (auto top = bigstruct.tops.begin(); top != bigstruct.tops.end(); ++top){
        totalpy += (*top)->calc.total_py();
    }
    return totalpy;
}

template<typename T> void pr1(T t){
    cout << left << setw(24) << setfill(' ') << t;
}

template<typename T> void pr(T t1, T t2, T t3, T t4){
    pr1(t1);
    pr1(t2);
    pr1(t3);
    pr1(t4);
    //pr1(t5);
}

inline void print_top_input(top_system & top){

    //Input values
    cout << endl;
    cout << "Input values:" << endl << endl;

    cout << "leptonic = " << top.input.leptonic << endl << endl;

    //Polar coords
    pr1(" ");
    pr("pt", "phi", "eta", "mass");
    cout << endl;

    pr1("b-jet");
    pr(top.input.b_pt, top.input.b_phi, top.input.b_eta, top.input.b_m);
    cout << endl;

    pr1("Wd1");
    pr(top.input.Wd1_pt, top.input.Wd1_phi, top.input.Wd1_eta, top.input.Wd1_m);
    cout << endl;

    pr1("Wd2");
    pr(top.input.Wd2_pt, top.input.Wd2_phi, top.input.Wd2_eta, top.input.Wd2_m);
    cout << endl;

    cout << endl;

    //Cartesian coords
    pr1(" ");
    pr("px", "py", "pz", "energy");
    cout << endl;

    pr1("b-jet");
    pr(top.input.b_px, top.input.b_py, top.input.b_pz, top.input.b_e);
    cout << endl;

    pr1("Wd1");
    pr(top.input.Wd1_px, top.input.Wd1_py, top.input.Wd1_pz, top.input.Wd1_e);
    cout << endl;

    pr1("Wd2");
    pr(top.input.Wd2_px, top.input.Wd2_py, top.input.Wd2_pz, top.input.Wd2_e);
    cout << endl;

    cout << endl;

/*    cout << " " << "\tpt       " << "\t\tphi     " << "\t\teta     " << "\t\tmass     " << endl;
    cout << "b-jet" << "\t" << top.input.b_pt << "\t\t" << top.input.b_phi
         << "\t\t" << top.input.b_eta << "\t\t" << top.input.b_m << endl;
*/
    return;

}

inline void print_top_current(top_system & top){

    //Input values
    cout << endl;
    cout << "Current values:" << endl << endl;

    cout << "leptonic = " << top.input.leptonic << endl << endl;

    //Polar coords
    pr1(" ");
    pr("pt", "phi", "eta", "input mass");
    cout << endl;

    pr1("b-jet");
    pr(top.calc.b_pt(), top.calc.b_phi(), top.calc.b_eta(), top.input.b_m);
    cout << endl;

    pr1("Wd1");
    pr(top.calc.Wd1_pt(), top.calc.Wd1_phi(), top.calc.Wd1_eta(), top.input.Wd1_m);
    cout << endl;

    pr1("Wd2");
    pr(top.calc.Wd2_pt(), top.calc.Wd2_phi(), top.calc.Wd2_eta(), top.input.Wd2_m);
    cout << endl;

    cout << endl;

    //Cartesian coords
    pr1(" ");
    pr("px", "py", "pz", "energy");
    cout << endl;

    pr1("b-jet");
    pr(top.calc.b_px(), top.calc.b_py(), top.calc.b_pz(), top.calc.b_e());
    cout << endl;

    pr1("Wd1");
    pr(top.calc.Wd1_px(), top.calc.Wd1_py(), top.calc.Wd1_pz(), top.calc.Wd1_e());
    cout << endl;

    pr1("Wd2");
    pr(top.vars.Wd2_px, top.vars.Wd2_py, top.vars.Wd2_pz, top.calc.Wd2_e());
    cout << endl;

    cout << endl;

    //Masses of parents
    cout << "mW = " << top.calc.mW() << endl;
    cout << "mTop = " << top.calc.mTop() << endl;
    cout << endl;

/*    cout << " " << "\tpt       " << "\t\tphi     " << "\t\teta     " << "\t\tmass     " << endl;
    cout << "b-jet" << "\t" << top.input.b_pt << "\t\t" << top.input.b_phi
         << "\t\t" << top.input.b_eta << "\t\t" << top.input.b_m << endl;
*/
    return;

}

inline void print_top_chi2s_current(top_system & top){

    //Input values
    cout << endl;
    cout << "Current chi2 values:" << endl << endl;

    cout << "leptonic = " << top.input.leptonic << endl << endl;

    cout << "b chi2 = " << top.calc.b_chi2() << endl;
    cout << "Wd1 chi2 = " << top.calc.Wd1_chi2() << endl;
    cout << "Wd2 chi2 = " << top.calc.Wd2_chi2() << endl;
    cout << "mW chi2 = " << top.calc.mW_chi2() << endl;
    cout << "mTop chi2 = " << top.calc.mTop_chi2() << endl;

    cout << endl;

    return;

}

inline void print_top_all(top_system &top){
    cout << setprecision(10) << " " << endl;
    print_top_input(top);
    print_top_current(top);
    print_top_chi2s_current(top);

    return;
}


}//END namespace commonstruct

#endif
