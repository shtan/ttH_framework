#ifndef maps_h
#define maps_h

#include <fstream>
#include <vector>
#include <map>
#include <TH1D.h>
#include <TCanvas.h>

using namespace std;

//namespace maps
//{

/*vector<string> particles = {"Bottom_1", "Wd11", "Wd12", "W1", "Top_1",
                            "Bottom_2", "Wd21", "Wd22", "W2", "Top_2",
                            "b1_from_H", "b2_from_H", "Higgs"};

vector<string> variables = {"Pt", "Phi", "Eta", "M"};*/

extern vector<string> particles;
extern vector<string> variables;
extern vector<string> fitstatus;
extern vector<string> dataset;
extern vector<string> datasetdiff;
extern vector<string> singleints;
extern vector<string> singledoubles;
extern vector<string> chi2s;
extern vector<string> diffvals;
extern vector<string> variablesCar;
extern vector<string> variablesH;
extern vector<string> variablesT;

typedef map<string, double> dmap1;
typedef map<string, dmap1> dmap2;
typedef map<string, dmap2> dmap3;

typedef map<string, int> imap1;

typedef map<string, ofstream*> fmap1;
typedef map<string, fmap1> fmap2;
typedef map<string, fmap2> fmap3;
typedef map<string, fmap3> fmap4;

typedef map<string, vector<double> > vdmap1;
typedef map<string, vdmap1> vdmap2;
typedef map<string, vdmap2> vdmap3;
typedef map<string, vdmap3> vdmap4;

typedef map<string, vector<int> > vimap1;
typedef map<string, vimap1> vimap2;

typedef map<string, TH1D *> hmap1;
typedef map<string, hmap1> hmap2;
typedef map<string, hmap2> hmap3;
typedef map<string, hmap3> hmap4;

typedef map<string, TCanvas*> cmap1;
typedef map<string, cmap1> cmap2;
typedef map<string, cmap2> cmap3;

//}

#endif
