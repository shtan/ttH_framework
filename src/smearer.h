#include <random>
#include <vector>

using namespace std;

struct resolution {
    double pT;  // relative resolution
    double phi; // resolution in rad.
    double eta; // resolution in rad.
};

class smearer
{
  public:
    smearer();

    std::default_random_engine generator;
    std::normal_distribution<double> normal;

    resolution j_res;
    resolution l_res;

    int j_max_id;

    void smear(const vector<pair<int, double *>> &,
               vector<pair<int, double *>> &);

  private:
    inline void smear_(const resolution &res, const double *p, const double m,
                       double *);
    inline void smear_givemass_(const resolution &res, const double *p, const double m,
                       double *);
    inline void resize(vector<pair<int, double *>> &p, const unsigned int sz);
};
