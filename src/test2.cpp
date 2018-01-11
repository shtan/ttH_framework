#ifndef test2_cpp
#define test2_cpp

#include <iostream>
#include <fstream>
#include <vector>
#include "LHEF_lite.h"
#include "analyzers_2.h"
#include "maps.h"
#include "scatter_plotter.h"

#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPDF.h>

using namespace std;
//using namespace maps;

void initialise_vecs();
void write_vecs(string);
void open_files_new(string);
void close_files_new();
void open_files(fmap2&, string, string);
void close_files(fmap2&);
void open_file(ofstream&, string, string);
void close_file(ofstream&);
void scatter_plot(vector<double>&, vector<double>&, string);
void scatter_plot_m(vector<double>&, vector<double>&, vector<double>&, string);
void read_files();
void plot();
vector<string> particles;
vector<string> variables;
vector<string> fitstatus;
vector<string> dataset;
vector<string> datasetdiff;
vector<string> singleints;
vector<string> singledoubles;
vector<string> chi2s;
vector<string> diffvals;

fmap4 file_diff_part_var;
fmap4 file_data_part_var;
fmap2 file_singleint;
fmap2 file_singledouble;
fmap2 file_chisquares;
fmap3 file_diff_diffvals;

vdmap4 vec_diff_part_var;
vdmap4 vec_data_part_var;
vimap2 vec_singleint;
vdmap2 vec_singledouble;
vdmap2 vec_chisquares;
vdmap3 vec_diff_diffvals;

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

    fitstatus = {"all", "converged", "failed"};

    dataset = {"gen", "smeared", "best"};

    datasetdiff = {"best_gen", "smeared_gen", "best_smeared"};

    singleints = {"inner_status", "outer_status", "event_number"};
    singledoubles = {"inner_edm", "outer_edm"};

    chi2s = {"Bottom_1", "Wd11", "Wd12", "mW1", "mTop_1",
                 "Bottom_2", "Wd21", "Wd22", "mW2", "mTop_2",
                 "nontops"};

    diffvals = {"diff_chi2_measurable", "diff_chi2_masses", "diff_chi2_neutrino",
                "diff_chi2_measurableplusmasses", "diff_chi2_total"};

    LHEF_lite reader("ttbb_h_bbbbdue.lhe");
    //LHEF_lite reader("ttbb_h_bbbbdue_10.lhe");

    unsigned long neve = atoi(argv[1]);
    unsigned long max_ev = atoi(argv[2]) + 1;

    string path_suffix = to_string(neve) + "_" + to_string(max_ev - 1);

    vector<pair<int, double *>> p; // particles
    vector<pair<int, int>> moth_ID;


    /*fmap4 file_diff_part_var;
    fmap4 file_data_part_var;
    fmap2 file_singleint;
    fmap2 file_singledouble;
    fmap2 file_chisquares;
    fmap3 file_diff_diffvals;*/

    /*cout << "before blah" << endl;
        ROOT::Math::Minimizer *blah = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Minimize");
        cout << "after blah"<< endl;*/

    //open_files_new(path_suffix);
    initialise_vecs();

//    read_files();

//    plot();

    //scatter_plot(vec_singledouble["converged"]["outer_edm"], vec_diff_diffvals["converged"]["best_gen"]["diff_chi2_total"], "./scatterplots/20170220/diffchi2_total_edm_converged_test.pdf");

    //scatter_plot_m(vec_singledouble["converged"]["outer_edm"], vec_diff_diffvals["converged"]["best_gen"]["diff_chi2_total"], vec_diff_diffvals["converged"]["smeared_gen"]["diff_chi2_total"], "./scatterplots/20170220/diffchi2_total_minus_edm_converged.pdf");

    //scatter_plot_m(vec_singledouble["converged"]["outer_edm"], vec_diff_part_var["converged"]["best_gen"]["Top_1"]["M"], vec_diff_part_var["converged"]["smeared_gen"]["Top_1"]["M"], "./scatterplots/20170220/diffchi2_total_minus_edm_converged.pdf");

    /*cout << "before blah" << endl;
        ROOT::Math::Minimizer *blah = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Minimize");
        cout << "after blah"<< endl;
    exit(0);*/
    //close_files_new();
    //exit(0);


    /*//fmap2 outfiles_best_gen;
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
    ofstream outfile_best_smeared_diffchi2_failed;*/

    /*//open_files(outfiles_best_gen, "best_gen");
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
    //outfile.open ( "./results/fit_output.txt" );*/

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

        a.analz(p, moth_ID, neve, vec_diff_part_var, vec_data_part_var, vec_singleint,
                vec_singledouble, vec_chisquares, vec_diff_diffvals);

        /*a.analz(p, moth_ID, neve, file_diff_part_var, file_data_part_var, file_singleint,
                file_singledouble, file_chisquares, file_diff_diffvals);*/

        
        /*a.analz(p, moth_ID, neve, outfiles_best_gen_all, outfiles_smeared_gen_all,
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
                );*/
        //cout << "exited analz" << endl;
    ++neve;
    }

    //close_files_new();
    write_vecs(path_suffix);

    /*cout << "before close files" << endl;

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

    //outfile.close();*/

    for (auto i : p)
        delete i.second;

    return 0;
}

void initialise_vecs()
{

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (auto v = variables.begin(); v != variables.end(); ++v){
                    const string var = *v;
                    //vec_diff_part_var[fits][diff][part][var] = new vector<double>;
                    //vec_diff_part_var[fits][diff][part][var] = ;
                    (vec_diff_part_var[fits][diff][part][var]).clear();
                }
            }
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = dataset.begin(); d != dataset.end(); ++d){
            const string data = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (auto v = variables.begin(); v != variables.end(); ++v){
                    const string var = *v;
                    //vec_data_part_var[fits][data][part][var] = new vector<double>;
                    (vec_data_part_var[fits][data][part][var]).clear();
                }
            }
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto s = singleints.begin(); s != singleints.end(); ++s){
            const string single = *s;
            //vec_singleint[fits][single] = new vector<double>;
            vec_singleint[fits][single].clear();
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto s = singledoubles.begin(); s != singledoubles.end(); ++s){
            const string single = *s;
            //vec_singledouble[fits][single] = new vector<double>;
            vec_singledouble[fits][single].clear();
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto c = chi2s.begin(); c != chi2s.end(); ++c){
            const string chi = *c;
            //vec_chisquares[fits][chi] = new vector<double>;
            vec_chisquares[fits][chi].clear();
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;
            for (auto v = diffvals.begin(); v != diffvals.end(); ++v){
                const string val = *v;
                //vec_diff_diffvals[fits][diff][val] = new vector<double>;
                vec_diff_diffvals[fits][diff][val].clear();
            }
        }
    }

}

void write_vecs( string path_suffix)
{
    //string bigpath = path + path_suffix + "/";
    string bigpath = "/afs/cern.ch/work/s/shtan/private/topreco_20161213/20170703/test/" + path_suffix + "/";
    //outpath = "/afs/cern.ch/user/s/shtan/teststorage/";

    //fmap outfiles;

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (auto v = variables.begin(); v != variables.end(); ++v){
                    const string var = *v;
                    ofstream* tempfile = new ofstream;
                    file_diff_part_var[fits][diff][part][var] = tempfile;
                    file_diff_part_var[fits][diff][part][var]->open (
                            (bigpath + diff + "_" + part + "_" + var + "_" + fits + ".txt").c_str() );
                    cout << "vector size = " << (vec_diff_part_var[fits][diff][part][var]).size() << endl;
                    for (auto ite = (vec_diff_part_var[fits][diff][part][var]).begin(); 
                            ite != (vec_diff_part_var[fits][diff][part][var]).end(); ++ite){
                        *file_diff_part_var[fits][diff][part][var] << *ite << endl;
                    }
                    (file_diff_part_var[fits][diff][part][var])->close();
                    
                    cout << "wrote file " << fits << diff << part << var << endl;
                }
            }
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = dataset.begin(); d != dataset.end(); ++d){
            const string data = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (auto v = variables.begin(); v != variables.end(); ++v){
                    const string var = *v;
                    ofstream* tempfile = new ofstream;
                    file_data_part_var[fits][data][part][var] = tempfile;
                    file_data_part_var[fits][data][part][var]->open ( (bigpath + data + "_" + part + "_" + var + "_" + fits + ".txt").c_str() );
                    cout << "vector size = " << (vec_data_part_var[fits][data][part][var]).size() << endl;
                    for (auto ite = (vec_data_part_var[fits][data][part][var]).begin(); 
                            ite != (vec_data_part_var[fits][data][part][var]).end(); ++ite){
                        *file_data_part_var[fits][data][part][var] << *ite << endl;
                    }
                    (file_data_part_var[fits][data][part][var])->close();
                    
                    cout << "wrote file " << fits << data << part << var << endl;
                }
            }
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto s = singleints.begin(); s != singleints.end(); ++s){
            const string single = *s;
            ofstream * tempfile = new ofstream;
            file_singleint[fits][single] = tempfile;
            file_singleint[fits][single]->open ( (bigpath + single + "_" + fits + ".txt").c_str() );
            cout << "vector size = " << (vec_singleint[fits][single]).size() << endl;
            for (auto ite = (vec_singleint[fits][single]).begin(); 
                    ite != (vec_singleint[fits][single]).end(); ++ite){
                *file_singleint[fits][single] << *ite << endl;
            }
            (file_singleint[fits][single])->close();
            
            cout << "wrote file " << fits << single << endl;
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto s = singledoubles.begin(); s != singledoubles.end(); ++s){
            const string single = *s;
            ofstream *tempfile = new ofstream;
            file_singledouble[fits][single] = tempfile;
            file_singledouble[fits][single]->open ( (bigpath + single + "_" + fits + ".txt").c_str() );
            for (auto ite = (vec_singledouble[fits][single]).begin(); 
                    ite != (vec_singledouble[fits][single]).end(); ++ite){
                *file_singledouble[fits][single] << *ite << endl;
            }
            (file_singledouble[fits][single])->close();
            
            cout << "wrote file " << fits << single << endl;
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto c = chi2s.begin(); c != chi2s.end(); ++c){
            const string chi = *c;
            ofstream* tempfile = new ofstream;
            file_chisquares[fits][chi] = tempfile;
            file_chisquares[fits][chi]->open ( (bigpath + "chi2_" + chi + "_" + fits + ".txt").c_str() );
            for (auto ite = (vec_chisquares[fits][chi]).begin(); 
                    ite != (vec_chisquares[fits][chi]).end(); ++ite){
                *file_chisquares[fits][chi] << *ite << endl;
            }
            (file_chisquares[fits][chi])->close();
            
            cout << "wrote file " << fits << chi << endl;
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;
            for (auto v = diffvals.begin(); v != diffvals.end(); ++v){
                const string val = *v;
                ofstream* tempfile = new ofstream;
                file_diff_diffvals[fits][diff][val] = tempfile;
                file_diff_diffvals[fits][diff][val]->open ( (bigpath + diff + "_" + val + "_" + fits + ".txt").c_str() );
                cout << "opening file " << fits << diff << val << endl;
                for (auto ite = (vec_diff_diffvals[fits][diff][val]).begin(); 
                        ite != (vec_diff_diffvals[fits][diff][val]).end(); ++ite){
                    *file_diff_diffvals[fits][diff][val] << *ite << endl;
                }
                (file_diff_diffvals[fits][diff][val])->close();
                
                cout << "wrote file " << fits << diff << val << endl;
            }
        }
    }

}

void read_files()
{
    //string bigpath = path + path_suffix + "/";
    string bigpath = "/afs/cern.ch/work/s/shtan/private/topreco_20161213/20170213/new_file_system/0_499/";
    //outpath = "/afs/cern.ch/user/s/shtan/teststorage/";

    //fmap outfiles;

    string line;

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (auto v = variables.begin(); v != variables.end(); ++v){
                    const string var = *v;
                    ifstream infile((bigpath + diff + "_" + part + "_" + var + "_" + fits + ".txt").c_str() );
                    while (getline(infile, line)){
                        vec_diff_part_var[fits][diff][part][var].push_back(atof(line.c_str()));
                    }
                    infile.close();
                    
                    cout << "read file " << fits << diff << part << var << endl;
                }
            }
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = dataset.begin(); d != dataset.end(); ++d){
            const string data = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (auto v = variables.begin(); v != variables.end(); ++v){
                    const string var = *v;
                    ifstream infile((bigpath + data + "_" + part + "_" + var + "_" + fits + ".txt").c_str() );
                    while (getline(infile, line)){
                        vec_data_part_var[fits][data][part][var].push_back(atof(line.c_str()));
                    }
                    infile.close();
                    
                    cout << "read file " << fits << data << part << var << endl;
                }
            }
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto s = singleints.begin(); s != singleints.end(); ++s){
            const string single = *s;
            ifstream infile((bigpath + single + "_" + fits + ".txt").c_str() );
            while (getline(infile, line)){
                vec_singleint[fits][single].push_back(atoi(line.c_str()));
            }

            infile.close();
            
            cout << "read file " << fits << single << endl;
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto s = singledoubles.begin(); s != singledoubles.end(); ++s){
            const string single = *s;
            ifstream infile((bigpath + single + "_" + fits + ".txt").c_str() );
            while (getline(infile, line)){
                vec_singledouble[fits][single].push_back(atof(line.c_str()));
            }
            infile.close();
            
            cout << "read file " << fits << single << endl;
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto c = chi2s.begin(); c != chi2s.end(); ++c){
            const string chi = *c;
            ifstream infile((bigpath + "chi2_" + chi + "_" + fits + ".txt").c_str() );
            while (getline(infile, line)){
                vec_chisquares[fits][chi].push_back(atof(line.c_str()));
            }
            infile.close();
            
            cout << "read file " << fits << chi << endl;
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;
            for (auto v = diffvals.begin(); v != diffvals.end(); ++v){
                const string val = *v;
                ifstream infile((bigpath + diff + "_" + val + "_" + fits + ".txt").c_str() );
                while (getline(infile, line)){
                    vec_diff_diffvals[fits][diff][val].push_back(atof(line.c_str()));
                }
                infile.close();
                
                cout << "read file " << fits << diff << val << endl;
            }
        }
    }

}

void open_files_new( string path_suffix)
{
    //string bigpath = path + path_suffix + "/";
    string bigpath = "/afs/cern.ch/work/s/shtan/private/topreco_20161213/20170213/nonexisting_folder/" + path_suffix + "/";
    //outpath = "/afs/cern.ch/user/s/shtan/teststorage/";

    //fmap outfiles;

    cout << file_diff_part_var.max_size() << endl;
    int counter = 0;

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (auto v = variables.begin(); v != variables.end(); ++v){
                    const string var = *v;
                    ofstream* tempfile = new ofstream;
                    //tempfile.open (
                    //        (bigpath + diff + "_" + part + "_" + var + "_" + fits + ".txt").c_str() );

                    file_diff_part_var[fits][diff][part][var] = tempfile;
                    file_diff_part_var[fits][diff][part][var]->open (
                            (bigpath + diff + "_" + part + "_" + var + "_" + fits + ".txt").c_str() );
                    //    .open (
                    //        (bigpath + diff + "_" + part + "_" + var + "_" + fits + ".txt").c_str() );
                    *file_diff_part_var[fits][diff][part][var] << "ha" << endl;
            //cout <<  (outpath + prefix + "_" + part + "_" + var + ".txt").c_str() << endl;
                    cout << "wrote file " << fits << diff << part << var << endl;
                    counter++;
                }
            }
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = dataset.begin(); d != dataset.end(); ++d){
            const string data = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (auto v = variables.begin(); v != variables.end(); ++v){
                    const string var = *v;
                    ofstream* tempfile = new ofstream;
                    file_data_part_var[fits][data][part][var] = tempfile;
                    file_data_part_var[fits][data][part][var]->open ( (bigpath + data + "_" + part + "_" + var + "_" + fits + ".txt").c_str() );
                    *file_data_part_var[fits][data][part][var] << "ha" << endl;
                    cout << "wrote file " << fits << data << part << var << endl;
                    counter++;
                }
            }
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto s = singleints.begin(); s != singleints.end(); ++s){
            const string single = *s;
            ofstream * tempfile = new ofstream;
            file_singleint[fits][single] = tempfile;
            file_singleint[fits][single]->open ( (bigpath + single + "_" + fits + ".txt").c_str() );
            *file_singleint[fits][single] << "ha" << endl;
            cout << "wrote file " << fits << single << endl;
            counter++;
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto s = singledoubles.begin(); s != singledoubles.end(); ++s){
            const string single = *s;
            ofstream *tempfile = new ofstream;
            file_singledouble[fits][single] = tempfile;
            file_singledouble[fits][single]->open ( (bigpath + single + "_" + fits + ".txt").c_str() );
            *file_singledouble[fits][single] << "ha" << endl;
            cout << "wrote file " << fits << single << endl;
            counter++;
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto c = chi2s.begin(); c != chi2s.end(); ++c){
            const string chi = *c;
            ofstream* tempfile = new ofstream;
            file_chisquares[fits][chi] = tempfile;
            file_chisquares[fits][chi]->open ( (bigpath + "chi2_" + chi + "_" + fits + ".txt").c_str() );
            *file_chisquares[fits][chi] << "ha" << endl;
            cout << "wrote file " << fits << chi << endl;
            counter++;
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;
            for (auto v = diffvals.begin(); v != diffvals.end(); ++v){
                const string val = *v;
                ofstream* tempfile = new ofstream;
                file_diff_diffvals[fits][diff][val] = tempfile;
                file_diff_diffvals[fits][diff][val]->open ( (bigpath + diff + "_" + val + "_" + fits + ".txt").c_str() );
                cout << "opening file " << fits << diff << val << endl;
                *file_diff_diffvals[fits][diff][val] << "ha" << endl;
                cout << "wrote file " << fits << diff << val << endl;
                counter++;
            }
        }
    }
    cout << "counter = " << counter << endl;

}

void close_files_new()
{

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (auto v = variables.begin(); v != variables.end(); ++v){
                    const string var = *v; 
                    file_diff_part_var[fits][diff][part][var]->close();
                }
            }
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = dataset.begin(); d != dataset.end(); ++d){
            const string data = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (auto v = variables.begin(); v != variables.end(); ++v){
                    const string var = *v;
                    file_data_part_var[fits][data][part][var]->close();
                }
            }
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto s = singleints.begin(); s != singleints.end(); ++s){
            const string single = *s;
            file_singleint[fits][single]->close();
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto s = singledoubles.begin(); s != singledoubles.end(); ++s){
            const string single = *s;
            file_singledouble[fits][single]->close();
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto c = chi2s.begin(); c != chi2s.end(); ++c){
            const string chi = *c;
            file_chisquares[fits][chi]->close();
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;
            for (auto v = diffvals.begin(); v != diffvals.end(); ++v){
                const string val = *v;
                file_diff_diffvals[fits][diff][val]->close();
            }
        }
    }

}

void scatter_plot(vector<double> &vecx, vector<double> &vecy, string filename)
{
    cout << "scattering plot" << endl;
    int n = vecx.size();
    double arrx[n];
    double arry[n];
    for (int i=0; i<n; i++){
        arrx[i] = log10(vecx.at(i)+1e-20);
        arry[i] = log10(vecy.at(i)+1e-20);
    }

    cout << "created arrs" << endl;
    TCanvas* canvas = new TCanvas("canvas");
    TGraph *g = new TGraph(n, arrx, arry);
    g->SetMarkerStyle(2);

    cout << "getting axes" << endl;
    g->Draw("ap");
    g->GetXaxis()->SetTitle("Log(outer_edm)");
    g->GetYaxis()->SetTitle("Log(diffchi2_total)");
    //g->GetXaxis()->SetRangeUser(0,1);
    //g->GetYaxis()->SetRangeUser(0,20000);
    g->SetTitle("Converged");
    cout << "before draw" << endl;
    //g->Draw("ap");
    //canvas->Update();

    cout << "after draw" << endl;
    canvas->Print(filename.c_str());
    cout << "after print" << endl;

}

void scatter_plot_m(vector<double> &vecx, vector<double> &vecy1, vector<double> &vecy2, string filename)
{
    cout << "scattering plot" << endl;
    int n = vecx.size();
    int m = vecy1.size();
    int m2 = vecy2.size();
    cout << n << " " << m << " " << m2 << endl;
    double arrx[n];
    double arry[n];
    int counter = 0;
    for (int i=0; i<n; i++){
        arrx[i] = log10(vecx.at(i)+1e-20);
        if (abs(vecy1.at(i)) > abs(vecy2.at(i))){
            arry[i] = log10(vecy1.at(i) - vecy2.at(i)+1e-20);
            cout << abs(vecy1.at(i)) - abs(vecy2.at(i)) << endl;
        } else {
            counter++;
            //cout << abs(vecy1.at(i)) - abs(vecy2.at(i)) << endl;
        }
         //   arry[i] = -log10(
        //arry[i] = log10( vecy1.at(i)+1e-20);
    }
    cout << "n = " << n << endl;
    cout << "counter = " << counter << endl;

    cout << "created arrs" << endl;
    TCanvas* canvas = new TCanvas("canvas");
    TGraph *g = new TGraph(n, arrx, arry);
    g->SetMarkerStyle(2);

    cout << "getting axes" << endl;
    g->Draw("ap");
    g->GetXaxis()->SetTitle("Log(outer_edm)");
    g->GetYaxis()->SetTitle("best_gen_diffchi2_total - smeared_gen_diffchi2_total");
    //g->GetXaxis()->SetRangeUser(0,1);
    //g->GetYaxis()->SetRangeUser(0,20000);
    g->SetTitle("Converged");
    //cout << "before draw" << endl;
    g->Draw("ap");
    canvas->Update();

    cout << "after draw" << endl;
    //canvas->Print(filename.c_str());
    canvas->Print("test.pdf");
    cout << "after print" << endl;

}

void plot()
{

    /*for (auto v = (vec_data_part_var["converged"]["gen"]["Wd22"]["Pt"]).begin(); v != (vec_data_part_var["converged"]["gen"]["Wd22"]["Pt"]).end();v++){
        cout << *v << endl;
    }*/


    string plotpath = "./scatterplots/20170220/goodbad/";
    string ypart, yvar;
/*
//    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
//        const string fits = *fs;
//        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
//            const string diff = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (auto v = variables.begin(); v != variables.end(); ++v){
                    const string var = *v;
                    plotter::plotter plot;
                        plot.good_bad(vec_singledouble["converged"]["outer_edm"], 
                                vec_diff_part_var["converged"]["best_gen"][part][var],
                                vec_diff_part_var["converged"]["smeared_gen"][part][var],
                                plotpath + part + "_" + var + ".pdf",
                                part + "_" + var + ", converged",
                                "Log10(outer_edm)");
                }
            }
//        }
//    }*/

    plotpath = "./scatterplots/20170221/features_data/";
    ypart = "Top_1";
    yvar = "M";
    for (auto p = particles.begin(); p != particles.end(); ++p){
        const string part = *p;
        for (auto v = variables.begin(); v != variables.end(); ++v){
            const string var = *v;
            plotter::plotter plot;
                plot.good_bad(vec_data_part_var["converged"]["gen"][part][var], 
                        vec_diff_part_var["converged"]["best_gen"][ypart][yvar],
                        vec_diff_part_var["converged"]["smeared_gen"][ypart][yvar],
                        plotpath + "Performance_" + ypart + "_" + yvar + "_vs_" + part + "_" + var + ".pdf",
                        ypart + "_" + yvar + ", converged",
                        part + "_" + var + "_gen");
        }
    }
 
    plotpath = "./scatterplots/20170221/features_data2/";
    ypart = "Wd12";
    yvar = "Eta";
    for (auto p = particles.begin(); p != particles.end(); ++p){
        const string part = *p;
        for (auto v = variables.begin(); v != variables.end(); ++v){
            const string var = *v;
            plotter::plotter plot;
                plot.good_bad(vec_data_part_var["converged"]["gen"][part][var], 
                        vec_diff_part_var["converged"]["best_gen"][ypart][yvar],
                        vec_diff_part_var["converged"]["smeared_gen"][ypart][yvar],
                        plotpath + "Performance_" + ypart + "_" + yvar + "_vs_" + part + "_" + var + ".pdf",
                        ypart + "_" + yvar + ", converged",
                        part + "_" + var + "_gen");
        }
    }
 
    plotpath = "./scatterplots/20170221/features_data_smeared/";
    ypart = "Top_1";
    yvar = "M";
    for (auto p = particles.begin(); p != particles.end(); ++p){
        const string part = *p;
        for (auto v = variables.begin(); v != variables.end(); ++v){
            const string var = *v;
            plotter::plotter plot;
                plot.good_bad(vec_data_part_var["converged"]["smeared"][part][var], 
                        vec_diff_part_var["converged"]["best_gen"][ypart][yvar],
                        vec_diff_part_var["converged"]["smeared_gen"][ypart][yvar],
                        plotpath + "Performance_" + ypart + "_" + yvar + "_vs_" + part + "_" + var + ".pdf",
                        ypart + "_" + yvar + ", converged",
                        part + "_" + var + "_smeared");
        }
    }
 
    plotpath = "./scatterplots/20170221/features_data_smeared2/";
    ypart = "Wd12";
    yvar = "Eta";
    for (auto p = particles.begin(); p != particles.end(); ++p){
        const string part = *p;
        for (auto v = variables.begin(); v != variables.end(); ++v){
            const string var = *v;
            plotter::plotter plot;
                plot.good_bad(vec_data_part_var["converged"]["smeared"][part][var], 
                        vec_diff_part_var["converged"]["best_gen"][ypart][yvar],
                        vec_diff_part_var["converged"]["smeared_gen"][ypart][yvar],
                        plotpath + "Performance_" + ypart + "_" + yvar + "_vs_" + part + "_" + var + ".pdf",
                        ypart + "_" + yvar + ", converged",
                        part + "_" + var + "_smeared");
        }
    }
 

}


/*void open_files( fmap2 &outfiles, string prefix, string path_suffix )
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
}*/

#endif
