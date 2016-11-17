#include "ttHReco/ttHReco.h"
#include "ttHReco/ttbar_reco_C/topEventMinimizer.h"
#include "ttHReco/ttbar_reco_C/commonstruct.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"

//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> LV;
typedef TLorentzVector LV;

namespace recoc
{

struct parameters_ttRC {
    double m_t;
    double m_W;

    double SD_m_t;
    double SD_m_W;
};

class ttRC
{
  public:
    // main method for evaluation
    output reco(const input &, const parameters_ttRC &);

    // constructor
    ttRC();

    void daughter_to_parents(const input_x<4> &, parents<4> &);
    void input_to_output(const input&, output&);
    void met_to_neutrino(const double &, const double &, double [4]);

  private:
    topEventMinimizer *rc;
};

void Print(const input &);
void Print(const output &);
void Print_diff(const output &, const input &);
void Print_diff(const input &, const input &);
}
