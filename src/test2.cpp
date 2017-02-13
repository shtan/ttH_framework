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

void open_files(fmap2&, string, string);
void close_files(fmap2&);
void open_file(ofstream&, string, string);
void close_file(ofstream&);
vector<string> particles;
vector<string> variables;

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << "STARTINGEVENT ENDINGEVENT" << std::endl;
        return 1;
    }
    particles = {"Bottom_1", "Wd11", "Wd12", "W1", "Top_1",
                 "Bottom_2", "Wd21", "Wd22", "W2", "Top_2",
                 "b1_from_H", "b2_from_H", "Higgs"};

    variables = {"Pt", "Phi", "Eta", "M"};

    LHEF_lite reader("ttbb_h_bbbbdue.lhe");
    unsigned long neve = atoi(argv[1]);
    unsigned long max_ev = atoi(argv[2]) + 1;

    string path_suffix = to_string(neve) + "_" + to_string(max_ev - 1);

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
    ofstream outfile_inner_edm_all;
    ofstream outfile_outer_edm_all;
    ofstream outfile_smeared_gen_diffchi2_all;
    ofstream outfile_best_gen_diffchi2_all;
    ofstream outfile_best_smeared_diffchi2_all;

    fmap2 outfiles_best_converged;
    fmap2 outfiles_smeared_converged;
    fmap2 outfiles_gen_converged;
    ofstream outfile_inner_status_converged;
    ofstream outfile_outer_status_converged;
    ofstream outfile_event_number_converged;
    ofstream outfile_inner_edm_converged;
    ofstream outfile_outer_edm_converged;
    ofstream outfile_smeared_gen_diffchi2_converged;
    ofstream outfile_best_gen_diffchi2_converged;
    ofstream outfile_best_smeared_diffchi2_converged;

    fmap2 outfiles_best_failed;
    fmap2 outfiles_smeared_failed;
    fmap2 outfiles_gen_failed;
    ofstream outfile_inner_status_failed;
    ofstream outfile_outer_status_failed;
    ofstream outfile_event_number_failed;
    ofstream outfile_inner_edm_failed;
    ofstream outfile_outer_edm_failed;
    ofstream outfile_smeared_gen_diffchi2_failed;
    ofstream outfile_best_gen_diffchi2_failed;
    ofstream outfile_best_smeared_diffchi2_failed;

    //open_files(outfiles_best_gen, "best_gen");
    //open_files(outfiles_smeared_gen, "smeared_gen");
    open_files(outfiles_best_gen_all, "best_gen_all", path_suffix);
    open_files(outfiles_smeared_gen_all, "smeared_gen_all", path_suffix);
    open_files(outfiles_best_gen_converged, "best_gen_converged", path_suffix);
    open_files(outfiles_smeared_gen_converged, "smeared_gen_converged", path_suffix);
    open_files(outfiles_best_gen_failed, "best_gen_failed", path_suffix);
    open_files(outfiles_smeared_gen_failed, "smeared_gen_failed", path_suffix);

    open_files(outfiles_best_all, "best_all", path_suffix);
    open_files(outfiles_smeared_all, "smeared_all", path_suffix);
    open_files(outfiles_gen_all, "gen_all", path_suffix);
    open_file(outfile_inner_status_all, "inner_status_all", path_suffix);
    open_file(outfile_outer_status_all, "outer_status_all", path_suffix);
    open_file(outfile_event_number_all, "event_number_all", path_suffix);
    open_file(outfile_inner_edm_all, "inner_edm_all", path_suffix);
    open_file(outfile_outer_edm_all, "outer_edm_all", path_suffix);
    open_file(outfile_smeared_gen_diffchi2_all, "smeared_gen_diffchi2_all", path_suffix);
    open_file(outfile_best_gen_diffchi2_all, "best_gen_diffchi2_all", path_suffix);
    open_file(outfile_best_smeared_diffchi2_all, "best_smeared_diffchi2_all", path_suffix);

    open_files(outfiles_best_converged, "best_converged", path_suffix);
    open_files(outfiles_smeared_converged, "smeared_converged", path_suffix);
    open_files(outfiles_gen_converged, "gen_converged", path_suffix);
    open_file(outfile_inner_status_converged, "inner_status_converged", path_suffix);
    open_file(outfile_outer_status_converged, "outer_status_converged", path_suffix);
    open_file(outfile_event_number_converged, "event_number_converged", path_suffix);
    open_file(outfile_inner_edm_converged, "inner_edm_converged", path_suffix);
    open_file(outfile_outer_edm_converged, "outer_edm_converged", path_suffix);
    open_file(outfile_smeared_gen_diffchi2_converged, "smeared_gen_diffchi2_converged", path_suffix);
    open_file(outfile_best_gen_diffchi2_converged, "best_gen_diffchi2_converged", path_suffix);
    open_file(outfile_best_smeared_diffchi2_converged, "best_smeared_diffchi2_converged", path_suffix);

    open_files(outfiles_best_failed, "best_failed", path_suffix);
    open_files(outfiles_smeared_failed, "smeared_failed", path_suffix);
    open_files(outfiles_gen_failed, "gen_failed", path_suffix);
    open_file(outfile_inner_status_failed, "inner_status_failed", path_suffix);
    open_file(outfile_outer_status_failed, "outer_status_failed", path_suffix);
    open_file(outfile_event_number_failed, "event_number_failed", path_suffix);
    open_file(outfile_inner_edm_failed, "inner_edm_failed", path_suffix);
    open_file(outfile_outer_edm_failed, "outer_edm_failed", path_suffix);
    open_file(outfile_smeared_gen_diffchi2_failed, "smeared_gen_diffchi2_failed", path_suffix);
    open_file(outfile_best_gen_diffchi2_failed, "best_gen_diffchi2_failed", path_suffix);
    open_file(outfile_best_smeared_diffchi2_failed, "best_smeared_diffchi2_failed", path_suffix);

    //ofstream outfile;
    //outfile.open ( "./results/fit_output.txt" );

    ex1::analyzer2 a;
   for (unsigned long int0 = 0; int0 < neve; int0++){
        reader.Get_event_GM(p, moth_ID);
        a.smear_only(p);
        cout << "Skipping event " << int0 << endl;
    }

    //ex1::analyzer2 a;
    while (reader.Get_event_GM(p, moth_ID) && neve < max_ev) {
        cout << "ANALYSING EVENT " << neve << endl;
        //++neve;

        a.analz(p, moth_ID, neve, outfiles_best_gen_all, outfiles_smeared_gen_all,
                outfiles_best_gen_converged, outfiles_smeared_gen_converged,
                outfiles_best_gen_failed, outfiles_smeared_gen_failed,
                outfiles_best_all, outfiles_smeared_all, outfiles_gen_all,
                outfile_inner_status_all, outfile_outer_status_all, outfile_event_number_all,
                outfile_inner_edm_all, outfile_outer_edm_all,
                outfiles_best_converged, outfiles_smeared_converged, outfiles_gen_converged,
                outfile_inner_status_converged, outfile_outer_status_converged, outfile_event_number_converged,
                outfile_inner_edm_converged, outfile_outer_edm_converged,
                outfiles_best_failed, outfiles_smeared_failed, outfiles_gen_failed,
                outfile_inner_status_failed, outfile_outer_status_failed, outfile_event_number_failed,
                outfile_inner_edm_failed, outfile_outer_edm_failed,
                outfile_smeared_gen_diffchi2_all, outfile_best_gen_diffchi2_all, outfile_best_smeared_diffchi2_all,
                outfile_smeared_gen_diffchi2_converged, outfile_best_gen_diffchi2_converged, outfile_best_smeared_diffchi2_converged,
                outfile_smeared_gen_diffchi2_failed, outfile_best_gen_diffchi2_failed, outfile_best_smeared_diffchi2_failed
                );
        cout << "exited analz" << endl;
    ++neve;
    }

    cout << "before close files" << endl;

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
    close_file(outfile_inner_edm_all);
    close_file(outfile_outer_edm_all);
    cout << "before close first new file" << endl;
    close_file(outfile_smeared_gen_diffchi2_all);
    cout << "after close first new file" << endl;
    close_file(outfile_best_gen_diffchi2_all);
    close_file(outfile_best_smeared_diffchi2_all);

    close_files(outfiles_best_converged);
    close_files(outfiles_smeared_converged);
    close_files(outfiles_gen_converged);
    close_file(outfile_inner_status_converged);
    close_file(outfile_outer_status_converged);
    close_file(outfile_event_number_converged);
    close_file(outfile_inner_edm_converged);
    close_file(outfile_outer_edm_converged);
    close_file(outfile_smeared_gen_diffchi2_converged);
    close_file(outfile_best_gen_diffchi2_converged);
    close_file(outfile_best_smeared_diffchi2_converged);

    close_files(outfiles_best_failed);
    close_files(outfiles_smeared_failed);
    close_files(outfiles_gen_failed);
    close_file(outfile_inner_status_failed);
    close_file(outfile_outer_status_failed);
    close_file(outfile_event_number_failed);
    close_file(outfile_inner_edm_failed);
    close_file(outfile_outer_edm_failed);
    close_file(outfile_smeared_gen_diffchi2_failed);
    close_file(outfile_best_gen_diffchi2_failed);
    close_file(outfile_best_smeared_diffchi2_failed);

    //outfile.close();

    for (auto i : p)
        delete i.second;

    return 0;
}

void open_files( fmap2 &outfiles, string prefix, string path_suffix )
{
    string outpath = "/afs/cern.ch/work/s/shtan/private/topreco_20161213/20170207/diffchi2/" + path_suffix + "/";
    //outpath = "/afs/cern.ch/user/s/shtan/teststorage/";

    //fmap outfiles;

    for (auto p = particles.begin(); p != particles.end(); ++p){
        const string part = *p;
        for (auto v = variables.begin(); v != variables.end(); ++v){
            const string var = *v;
            
            outfiles[part][var].open ( (outpath + prefix + "_" + part + "_" + var + ".txt").c_str() );
            //cout <<  (outpath + prefix + "_" + part + "_" + var + ".txt").c_str() << endl;
        }
    }
}

void open_file( ofstream &file, string name, string path_suffix )
{
    string outpath = "/afs/cern.ch/work/s/shtan/private/topreco_20161213/20170207/diffchi2/" + path_suffix + "/";
//afs/cern.ch/work/s/shtan/private/topreco_20161213/20161213/56_59
    
    //outpath = "./results/test/thisfolder/";
    //outpath = "/afs/cern.ch/user/s/shtan/teststorage/";
    //string outpath = "./results/test/";
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
