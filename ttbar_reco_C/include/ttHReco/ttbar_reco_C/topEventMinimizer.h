#ifndef TOPEVENTMINIMIZER
#define TOPEVENTMINIMIZER

//#include "neutrinoSolutions.h"
#include "topSystemChiSquare.h"
//#include "hadronicTopSystemChiSquare.h"
//#include "leptonicTopSystemChiSquare.h"
#include "lightJetChiSquareMinimumSolver.h"
#include "commonstruct.h"
#include "TLorentzVector.h"
#include <map>

//using namespace commonstruct;

class topEventMinimizer
{
  private:
    commonstruct::big_struct & bigstruct;
    commonstruct::nontop_system & nontops;

    vector< topSystemChiSquare > topSysChiSqs_;


    lightJetChiSquareMinimumSolver nonTopChiSquare_;


    double maxConsideredChiSquareRoot_;

    ROOT::Math::Minimizer *innerMin_;
    ROOT::Math::Minimizer *outerMin_;

    //DON'T DELETE THIS
    //bool checkInputSizes();

    // void setBJets();
    //void setNonTopObjectCollections();
    // void setWDaughters();

    //  void initializeDeltas();
    //void initializeChiSquares();
    void Initialize_minimizers(ROOT::Math::Minimizer *&outer,
                               ROOT::Math::Minimizer *&inner);
    void initialize_best_outer_chiSquares();
    void reset_best_inner_chiSquares();

    void calcWDaughterEllipses();

    void setBestValues();

    int& debug_verbosity;

  public:

    //commonstruct::big_struct & bigstruct;


    // flags
    const bool debug = true;

    // status observables
    //int innerMinStatus;
    //int outerMinStatus;

    //double outerMin_Edm;

    // functions
    topEventMinimizer(commonstruct::big_struct&, int&);
    

    ~topEventMinimizer();


    void create_tops();

    void printTopConstituents();
    void printNonTopObjects();
    void checkRecoil_after_fit();

    double outerMinimizationOperator(const double *);
    double innerMinimizationOperator(const double *);

    void findStartingValues(int);
    void minimizeNonTopChiSquare();
    void minimizeTotalChiSquare();

    TLorentzVector get_b(int);
    TLorentzVector get_Wd1(int);
    TLorentzVector get_Wd2(int);
    TLorentzVector get_W(int);
    TLorentzVector get_top(int);
    TLorentzVector get_nontop_object(int);
    TLorentzVector get_higgs();
    double get_best_mW_chi2(int);
    double get_best_mTop_chi2(int);
    double get_best_b_chi2(int);
    double get_best_Wd1_chi2(int);
    double get_best_Wd2_chi2(int);
    double get_best_nontop_chi2();
    double get_best_total_had_chi2();
    double get_best_total_mTop_chi2();
    double get_best_total_topsys_chi2();
    double get_best_total_chi2();
    int get_inner_min_status();
    int get_outer_min_status();
    double get_inner_edm();
    double get_outer_edm();
    void plot_chi2s();
    void plot_chi2_outer(int, double&, const double&, const double&, double&, string);
    void plot_Wd22_chi2();



};

#endif
