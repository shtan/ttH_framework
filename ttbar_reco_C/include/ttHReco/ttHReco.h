#include <vector>

// reco core
namespace recoc
{

// struct parameters {
//     bool debug;
//     diagnostics diag; // if debug = true
//
//     parameters_CU pars1;
//     parameters_KIT pars2;
//     parameters_RWTH pars3;
//
//     parameters() : debug(false);
// };
//
// struct parameters_CU {
// };
//
// struct parameters_KIT {
// };
//
// struct parameters_RWTH {
// };
//
// class diagnostics
// {
// };

// polar: pT, phi, eta, E
template <unsigned int sz> struct input_x {
    //     const unsigned int sz;
    // top 1
    bool t1_lep; // leptonic
    double b1[sz];
    double d11[sz]; // possibly leptop
    double d12[sz]; // unused if neutrino

    // top 2
    bool t2_lep; // leptonic
    double b2[sz];
    double d21[sz]; // possibly lepton
    double d22[sz]; // unused if neutrino

    // bs from H
    double bH1[sz];
    double bH2[sz];

    std::vector<double *> p_others;

    //double MET_px;
    //double MET_py;
};

template <unsigned int sz> struct neutrino {
    //for generated neutrino storage
    double n1[sz];
    double n2[sz];
};

struct input {
    input_x<4> p;
    input_x<3> SD;
    double MET_px;
    double MET_py;
    neutrino<4> neu_p;
};

template <unsigned int sz> struct parents {
    //top 1
    double t1[sz];
    double w1[sz];

    //top 2
    double t2[sz];
    double w2[sz];

    //Higgs
    double h[sz];
};

struct chi_squareds {
    //top1
    double mt1;
    double mw1;
    double b1;
    double d11;
    double d12;

    //top 2
    double mt2;
    double mw2;
    double b2;
    double d21;
    double d22;

    //nontops
    double nontop_objs;
};

struct output {
    input_x<4> p;
    //input_x<1> particle_chi2s;
    parents<4> parents_p;
    //parents<1> parents_m_chi2s;
    chi_squareds chi2s;

    int inner_min_status;
    int outer_min_status;
    double weight;
};

class ttH_reconstructor
{
  public:
    // settings
    //     parameters params;

    // main method for evaluation
    output reco(const input &);

    // constructor
    ttH_reconstructor();

  private:
    // otherwise, we could merge inside
    //     algo_CU a1;
    //     algo_KIT a2;
    //     algo_RWTH a3;
};
}
