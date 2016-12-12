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
void open_file(ofstream&, string);
void close_file(ofstream&);
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
    const unsigned long max_ev = 2;

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

    fmap2 outfiles_best_all;
    fmap2 outfiles_smeared_all;
    fmap2 outfiles_gen_all;
    ofstream outfile_inner_status_all;
    ofstream outfile_outer_status_all;
    ofstream outfile_event_number_all;

    fmap2 outfiles_best_converged;
    fmap2 outfiles_smeared_converged;
    fmap2 outfiles_gen_converged;
    ofstream outfile_inner_status_converged;
    ofstream outfile_outer_status_converged;
    ofstream outfile_event_number_converged;

    fmap2 outfiles_best_failed;
    fmap2 outfiles_smeared_failed;
    fmap2 outfiles_gen_failed;
    ofstream outfile_inner_status_failed;
    ofstream outfile_outer_status_failed;
    ofstream outfile_event_number_failed;

    //open_files(outfiles_best_gen, "best_gen");
    //open_files(outfiles_smeared_gen, "smeared_gen");
    open_files(outfiles_best_gen_all, "best_gen_all");
    open_files(outfiles_smeared_gen_all, "smeared_gen_all");
    open_files(outfiles_best_gen_converged, "best_gen_converged");
    open_files(outfiles_smeared_gen_converged, "smeared_gen_converged");
    open_files(outfiles_best_gen_failed, "best_gen_failed");
    open_files(outfiles_smeared_gen_failed, "smeared_gen_failed");

    open_files(outfiles_best_all, "best_all");
    open_files(outfiles_smeared_all, "smeared_all");
    open_files(outfiles_gen_all, "gen_all");
    open_file(outfile_inner_status_all, "inner_status_all");
    open_file(outfile_outer_status_all, "outer_status_all");
    open_file(outfile_event_number_all, "event_number_all");

    open_files(outfiles_best_converged, "best_converged");
    open_files(outfiles_smeared_converged, "smeared_converged");
    open_files(outfiles_gen_converged, "gen_converged");
    open_file(outfile_inner_status_converged, "inner_status_converged");
    open_file(outfile_outer_status_converged, "outer_status_converged");
    open_file(outfile_event_number_converged, "event_number_converged");

    open_files(outfiles_best_failed, "best_failed");
    open_files(outfiles_smeared_failed, "smeared_failed");
    open_files(outfiles_gen_failed, "gen_failed");
    open_file(outfile_inner_status_failed, "inner_status_failed");
    open_file(outfile_outer_status_failed, "outer_status_failed");
    open_file(outfile_event_number_failed, "event_number_failed");

    //ofstream outfile;
    //outfile.open ( "./results/fit_output.txt" );

    ex1::analyzer2 a;
    while (reader.Get_event_GM(p, moth_ID) && neve < max_ev) {
        cout << "ANALYSING EVENT " << neve << endl;
        ++neve;

        a.analz(p, moth_ID, neve, outfiles_best_gen_all, outfiles_smeared_gen_all,
                outfiles_best_gen_converged, outfiles_smeared_gen_converged,
                outfiles_best_gen_failed, outfiles_smeared_gen_failed,
                outfiles_best_all, outfiles_smeared_all, outfiles_gen_all,
                outfile_inner_status_all, outfile_outer_status_all, outfile_event_number_all,
                outfiles_best_converged, outfiles_smeared_converged, outfiles_gen_converged,
                outfile_inner_status_converged, outfile_outer_status_converged, outfile_event_number_converged,
                outfiles_best_failed, outfiles_smeared_failed, outfiles_gen_failed,
                outfile_inner_status_failed, outfile_outer_status_failed, outfile_event_number_failed
                );
    }

    close_files(outfiles_best_gen_all);
    close_files(outfiles_smeared_gen_all);
    close_files(outfiles_best_gen_converged);
    close_files(outfiles_smeared_gen_converged);
    close_files(outfiles_best_gen_failed);
    close_files(outfiles_smeared_gen_failed);

    close_files(outfiles_best_all);
    close_files(outfiles_smeared_all);
    close_files(outfiles_gen_all);
    close_file(outfile_inner_status_all);
    close_file(outfile_outer_status_all);
    close_file(outfile_event_number_all);

    close_files(outfiles_best_converged);
    close_files(outfiles_smeared_converged);
    close_files(outfiles_gen_converged);
    close_file(outfile_inner_status_converged);
    close_file(outfile_outer_status_converged);
    close_file(outfile_event_number_converged);

    close_files(outfiles_best_failed);
    close_files(outfiles_smeared_failed);
    close_files(outfiles_gen_failed);
    close_file(outfile_inner_status_failed);
    close_file(outfile_outer_status_failed);
    close_file(outfile_event_number_failed);

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

void open_file( ofstream &file, string name )
{
    string outpath = "./results/test/";
    file.open( (outpath + name + ".txt").c_str() );
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

void close_file( ofstream &file )
{
    file.close();
}

#endif
