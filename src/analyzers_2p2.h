#include "analyzers.h"
#include "analyzers_1.h"
// #include "smearer.h" // in "analyzers_1.h"
// #include "identifier.h"  // in "analyzers_1.h"
// #include "MEKD.h"    // in "analyzers_1.h"
#include "ttHReco/ttbar_reco_C.h"
#include <fstream>

using namespace std;

namespace ex1
{

class analyzer2p2
{
  public:
    analyzer2p2();
    ~analyzer2p2();

    void analz(const pvec &, const movec &, ofstream &, ofstream &, ofstream &);
    inline void eval_MEs_ttbb(const ttbarX &ttH, vector<double> &MEs);

    smearer smr;
    pvec ps; // smeared particles
    ex1::identifier_GM ident;

  private:
    mekd::MEKD *mes;
    mekd::input in;
    vector<int> event_id;
    vector<double *> event_p;
    vector<double> MEs;
    ttbarX temp;
    recoc::input generated;
    recoc::input in_2_RC;
    recoc::parameters_ttRC params;
    recoc::ttRC reco_C;
};
}
