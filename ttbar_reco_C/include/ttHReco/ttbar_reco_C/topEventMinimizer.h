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
    double get_best_total_had_chi2();
    double get_best_total_mTop_chi2();
    double get_best_total_topsys_chi2();
    double get_best_total_chi2();





};

#endif
