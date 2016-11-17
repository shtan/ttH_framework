#ifndef analyzers2_h
#define analyzers2_h

#include "analyzers.h"
#include "smearer.h"
#include "identifier.h"
#include "ttHReco/ttbar_reco_C.h"
#include "maps.h"
#include <fstream>

using namespace std;

namespace ex1
{

class analyzer2
{
  public:
    analyzer2();
    ~analyzer2();

    //void analz(const pvec &, const movec &, ofstream &);
    void analz(const pvec &, const movec &, fmap2 &, fmap2 &, fmap2 &, fmap2 &, fmap2 &, fmap2 &);
    void input_to_output(const recoc::input &, recoc::output &);
    void calc_diff(const recoc::output &, const recoc::output &, dmap2 &);
    void principal_angle(double&);
    void one_diff(const double [4], const double [4], dmap2 &, string);
    void write_diff(dmap2 &, fmap2 &);

    smearer smr;
    vector<pair<int, double *>> ps; // smeared particles
    ex1::identifier_GM ident;
    recoc::input generated;
    recoc::input in_2_RC;
    recoc::parameters_ttRC params;
    recoc::ttRC reco_C;

    dmap2 best_gen;
    dmap2 smeared_gen;
    /*dmap2 best_gen_all;
    dmap2 smeared_gen_all;
    dmap2 best_gen_converged;
    dmap2 smeared_gen_converged;
    dmap2 best_gen_failed;
    dmap2 smeared_gen_failed;*/

/*    vector<string> particles = {"Bottom_1", "Wd11", "Wd12", "W1", "Top_1",
                                "Bottom_2", "Wd21", "Wd22", "W2", "Top_2",
                                "b1_from_H", "b2_from_H", "Higgs"};
    
    vector<string> variables = {"Pt", "Phi", "Eta", "M"};

    typedef map<string, double> dmap1;
    typedef map<string, dmap1> dmap2;

    dmap2 best_gen;
    dmap2 smeared_gen;*/

  private:
    void Setup_Maps();

};
}

#endif
