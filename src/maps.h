#ifndef maps_h
#define maps_h

#include <fstream>
#include <vector>
#include <map>

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
//double blah;

typedef map<string, double> dmap1;
typedef map<string, dmap1> dmap2;
typedef map<string, dmap2> dmap3;

typedef map<string, int> imap1;

//dmap2 best_gen;
//dmap2 smeared_gen;

typedef map<string, ofstream*> fmap1;
typedef map<string, fmap1> fmap2;
typedef map<string, fmap2> fmap3;
typedef map<string, fmap3> fmap4;


//}

#endif