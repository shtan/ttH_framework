#ifndef test2_cpp
#define test2_cpp

#include <iostream>
#include <fstream>
#include <vector>
#include "LHEF_lite.h"
#include "analyzers_2.h"
#include "maps.h"

using namespace std;
//using namespace maps;

void open_files(fmap2&, string);
void close_files(fmap2&);
vector<string> particles;
vector<string> variables;

int main()
{
    particles = {"Bottom_1", "Wd11", "Wd12", "W1", "Top_1",
                 "Bottom_2", "Wd21", "Wd22", "W2", "Top_2",
                 "b1_from_H", "b2_from_H", "Higgs"};

    variables = {"Pt", "Phi", "Eta", "M"};

    LHEF_lite reader("ttbb_h_bbbbdue.lhe");
    unsigned long neve = 0;
    const unsigned long max_ev = 10;

    vector<pair<int, double *>> p; // particles
    vector<pair<int, int>> moth_ID;

    //fmap2 outfiles_best_gen;
    //fmap2 outfiles_smeared_gen;
    fmap2 outfiles_best_gen_all;
    fmap2 outfiles_smeared_gen_all;
    fmap2 outfiles_best_gen_converged;
    fmap2 outfiles_smeared_gen_converged;
    fmap2 outfiles_best_gen_failed;
    fmap2 outfiles_smeared_gen_failed;
    //open_files(outfiles_best_gen, "best_gen");
    //open_files(outfiles_smeared_gen, "smeared_gen");
    open_files(outfiles_best_gen_all, "best_gen_all");
    open_files(outfiles_smeared_gen_all, "smeared_gen_all");
    open_files(outfiles_best_gen_converged, "best_gen_converged");
    open_files(outfiles_smeared_gen_converged, "smeared_gen_converged");
    open_files(outfiles_best_gen_failed, "best_gen_failed");
    open_files(outfiles_smeared_gen_failed, "smeared_gen_failed");

    //ofstream outfile;
    //outfile.open ( "./results/fit_output.txt" );

    ex1::analyzer2 a;
    while (reader.Get_event_GM(p, moth_ID) && neve < max_ev) {
        cout << "ANALYSING EVENT " << neve << endl;
        ++neve;

        a.analz(p, moth_ID, outfiles_best_gen_all, outfiles_smeared_gen_all,
                outfiles_best_gen_converged, outfiles_smeared_gen_converged,
                outfiles_best_gen_failed, outfiles_smeared_gen_failed);
    }

    close_files(outfiles_best_gen_all);
    close_files(outfiles_smeared_gen_all);
    close_files(outfiles_best_gen_converged);
    close_files(outfiles_smeared_gen_converged);
    close_files(outfiles_best_gen_failed);
    close_files(outfiles_smeared_gen_failed);

    //outfile.close();

    for (auto i : p)
        delete i.second;

    return 0;
}

void open_files( fmap2 &outfiles, string prefix )
{
    string outpath = "./results/test/";

    //fmap outfiles;

    for (auto p = particles.begin(); p != particles.end(); ++p){
        const string part = *p;
        for (auto v = variables.begin(); v != variables.end(); ++v){
            const string var = *v;
            
            outfiles[part][var].open ( (outpath + prefix + "_" + part + "_" + var + ".txt").c_str() );
        }
    }
}

void close_files( fmap2 &outfiles )
{
    for (auto p = particles.begin(); p != particles.end(); ++p){
        const string part = *p;
        for (auto v = variables.begin(); v != variables.end(); ++v){
            const string var = *v;
            
            outfiles[part][var].close();
        }
    }
}

#endif
