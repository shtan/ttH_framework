#include "analyzers.h"
#include "smearer.h"
#include "identifier.h"
#include "ttHReco/ttbar_reco_C.h"

using namespace std;

namespace ex1
{

class analyzer2
{
  public:
    analyzer2();
    ~analyzer2();

    void analz(const pvec &, const movec &);

    smearer smr;
    vector<pair<int, double *>> ps; // smeared particles
    ex1::identifier_GM ident;
    recoc::input in_2_RC;
    recoc::parameters_ttRC params;
    recoc::ttRC reco_C;
};
}