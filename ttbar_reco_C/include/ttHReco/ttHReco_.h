// reco core
namespace recoc
{

struct parameters {
    bool debug;
    diagnostics diag; // if debug = true

    parameters_CU pars1;
    parameters_KIT pars2;
    parameters_RWTH pars3;

    parameters() : debug(false);
};

struct parameters_CU {
};

struct parameters_KIT {
};

struct parameters_RWTH {
};

class diagnostics
{
};

struct input {
    input_x p;
    input_x SD;

    input() : p(4), SD(3);
}

struct input_x {
    const unsigned int sz;
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

    vector<double *> p_others;

    input_x(unsigned int D) : sz(D);
};

struct output {
    input_x p;
    double weight;

    output() : p(4);
};

class ttH_reconstructor
{
  public:
    // settings
    parameters params;

    // main method for evaluation
    output reco(const input &);

    // constructor
    ttH_reconstructor();

  private:
    // otherwise, we could merge inside
    algo_CU a1;
    algo_KIT a2;
    algo_RWTH a3;
};
}