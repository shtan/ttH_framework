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
    void smear_only(const pvec &);
    void analz(const pvec &, const movec &, long unsigned int &, vdmap4 &, vdmap4 &, vimap2 &, vdmap2 &, vdmap2 &, vdmap3 & );
    void input_to_output(const recoc::input &, recoc::output &);
    void calc_diff(const recoc::output &, const recoc::output &, dmap2 &);
    void principal_angle(double&);
    void one_diff(const double [4], const double [4], dmap2 &, string);
    void push_diff(dmap2 &, vdmap2 &);
    void write_diff(dmap2 &, fmap2 &);
    void write_big(string, recoc::output &, long unsigned int &, fmap4&, fmap4&, fmap2&, fmap2&, fmap2&, fmap3&);
    void push_big(string, recoc::output &, long unsigned int &, vdmap4&, vdmap4&, vimap2&, vdmap2&, vdmap2&, vdmap3&);
    void dmap2_converter(const recoc::output &, dmap2 &);
    void single_converter(const double [4], dmap2 &, string);
    void push_int(int &, vector<int>&);
    void push_int(long unsigned int &, vector<int>&);
    void push_int(double &, vector<double>&);
    void write_int(int &, ofstream *);
    void write_int(long unsigned int &, ofstream *);
    void write_int(double &, ofstream *);
    double one_p_chi2(dmap2 &, string, const double [3]);
    double one_m_chi2(dmap2 &, dmap2 &, string, const double);
    double breitWignerErr(const double &, const double &, const double &);
    double diff_chi2(dmap2 &, dmap2 &, recoc::output &, string);
    void get_chi2s(recoc::output &, dmap1 &);

    smearer smr;
    vector<pair<int, double *>> ps; // smeared particles
    ex1::identifier_GM ident;
    recoc::input generated;
    recoc::input in_2_RC;
    recoc::parameters_ttRC params;
    recoc::ttRC reco_C;

    dmap3 diff_part_var;
    dmap3 data_part_var;
    imap1 singleint;
    dmap1 singledouble;
    dmap1 chisquares;
    dmap2 diff_diffvals;


  private:
    void Setup_Maps();

};
}

#endif
