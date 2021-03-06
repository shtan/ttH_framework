#ifndef test2_cpp
#define test2_cpp

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "LHEF_lite.h"
#include "analyzers_2.h"
#include "maps.h"
#include "scatter_plotter.h"

#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPDF.h>
#include <TPaveStats.h>
#include <TFile.h>
#include <TGaxis.h>

using namespace std;
//using namespace maps;

void perform_fit(const string, const string, unsigned long, unsigned long);
void perform_plot(const string);
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
void read_files(const string);
void plot();
void initialise_hists();
void convert_vecs_to_car();
void fill_hists();
void declare_canvases();
void plot_hists(const string);
void moveStatsBox(TH1D*);

vector<string> particles;
vector<string> variables;
vector<string> fitstatus;
vector<string> dataset;
vector<string> datasetdiff;
vector<string> singleints;
vector<string> singledoubles;
vector<string> chi2s;
vector<string> diffvals;
vector<string> variablesCar;
vector<string> variablesH;
vector<string> variablesT;

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

vdmap4 vec_data_part_varcar; //For Cartesian variables and E

hmap4 hist_diff_part_var;
cmap3 canvs_part_var;

int main(int argc, char* argv[])
{
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " STARTINGEVENT ENDINGEVENT inputfile outputdir task" << std::endl;
        return 1;
    }
    
    //Set up name lists
    particles = {"Bottom_1", "Wd11", "Wd12", "W1", "Top_1",
                 "Bottom_2", "Wd21", "Wd22", "W2", "Top_2",
                 "b1_from_H", "b2_from_H", "Higgs"};

    variables = {"Pt", "Phi", "Eta", "M"};
    
    variablesCar = {"Px", "Py", "Pz", "E"};
    
    variablesT = {"Pt", "Phi", "Eta", "M", "Px", "Py", "Pz", "E"};

    //for plotting histograms
    variablesH = {"Pt", "Phi", "Eta", "M", //for residuals
                 "Px", "Py", "Pz", "E",   // for residuals
                 "Pt_", "Phi_", "Eta_", "M_", // for resolutions
                 "Px_", "Py_", "Pz_", "E_"};  // for resolutions
    
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

    
    //Read in arguments

    unsigned long neve = atoi(argv[1]); //starting event
    unsigned long max_ev = atoi(argv[2]) + 1; //ending event
    const string input(argv[3]); //input file (for fit or fitplot task) or input directory (for plot task)
    const string outdir(argv[4]); //output directory
    const string task(argv[5]);
    // task can be:
    // fit: do fit, output results in text files
    // plot: read in text files, produce plots
    // fitplot: do fit, output results in text files, and make plots
    
    string path_suffix = to_string(neve) + "_" + to_string(max_ev - 1);
    
    if (task == "fit"){
        string outputdir = outdir + "/";
        perform_fit(input, outputdir, neve, max_ev);
    } else if (task == "fitplot") {
        string outputdir = outdir + "/";
        perform_fit(input, outputdir, neve, max_ev);
        perform_plot(outputdir);
    } else if (task == "plot") {
        read_files(input + "/");
        perform_plot(outdir + "/");
    }
    

    return 0;
}

void perform_fit(const string input, const string outdir, unsigned long neve,
                 unsigned long max_ev)
{
    LHEF_lite reader(input);

    //Folder to put output files into, according to starting and ending events
    //string path_suffix = to_string(neve) + "_" + to_string(max_ev - 1);

    vector<pair<int, double *>> p; // particles
    vector<pair<int, int>> moth_ID;
    
    initialise_vecs();
    
    ex1::analyzer2 a;
    
    //Hack to get reader to start reading from starting event
    for (unsigned long int0 = 0; int0 < neve; int0++){
        reader.Get_event_GM(p, moth_ID);
        a.smear_only(p);
        cout << "Skipping event " << int0 << endl;
    }

    //Do actual fits for given event range
    while (reader.Get_event_GM(p, moth_ID) && neve < max_ev) {
        cout << "ANALYSING EVENT " << neve << endl;

        a.analz(p, moth_ID, neve, vec_diff_part_var, vec_data_part_var, vec_singleint,
                vec_singledouble, vec_chisquares, vec_diff_diffvals);
        
        ++neve;
    }
    
    write_vecs(outdir);
    
    for (auto i : p)
        delete i.second;
}

void perform_plot(const string outdir)
{
    TFile* outfileplot = new TFile( (outdir + "plots/output.root").c_str(), "RECREATE");
    outfileplot->cd();
    initialise_hists();
    convert_vecs_to_car();
    fill_hists();
    declare_canvases();
    plot_hists( outdir + "plots/" );
    outfileplot->Write();
    outfileplot->Close();
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
                    (vec_data_part_var[fits][data][part][var]).clear();
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
                for (auto v = variablesCar.begin(); v != variablesCar.end(); ++v){
                    const string var = *v;
                    (vec_data_part_varcar[fits][data][part][var]).clear();
                }
            }
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto s = singleints.begin(); s != singleints.end(); ++s){
            const string single = *s;
            vec_singleint[fits][single].clear();
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto s = singledoubles.begin(); s != singledoubles.end(); ++s){
            const string single = *s;
            vec_singledouble[fits][single].clear();
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto c = chi2s.begin(); c != chi2s.end(); ++c){
            const string chi = *c;
            vec_chisquares[fits][chi].clear();
        }
    }

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;
            for (auto v = diffvals.begin(); v != diffvals.end(); ++v){
                const string val = *v;
                vec_diff_diffvals[fits][diff][val].clear();
            }
        }
    }

}

//void write_vecs( string path_suffix)
void write_vecs( string bigpath)
{
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
                    //cout << "vector size = " << (vec_diff_part_var[fits][diff][part][var]).size() << endl;
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
                    //cout << "vector size = " << (vec_data_part_var[fits][data][part][var]).size() << endl;
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
            //cout << "vector size = " << (vec_singleint[fits][single]).size() << endl;
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
                //cout << "opening file " << fits << diff << val << endl;
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

void initialise_hists()
{
    double lbound = -1;
    double rbound = 1;

    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (auto v = variablesH.begin(); v != variablesH.end(); ++v){
                    const string var = *v;
                    
                    //Set upper and lower histogram bounds according to variable
                    //For residual plots
                    if (var == "Pt" || var == "Px" || var == "Py" || var == "Pz") {
                        lbound = -150;
                        rbound = 150;
                    } else if (var == "Eta") {
                        lbound = -5;
                        rbound = 5;
                    } else if (var == "Phi") {
                        lbound = -5;
                        rbound = 5;
                    } else if (var == "M" || var == "E") {
                        lbound = -150;
                        rbound = 150;
                    } else {
                        // If resolution plots, use -1 to 1
                        lbound = -1;
                        rbound = 1;
                    }
                    
 
                    string hname;
                    if (diff == "best_gen") hname = "fit-gen";
                    else if (diff == "smeared_gen") hname = "smeared-gen";
                    else if (diff == "smeared_best") hname = "smeared-fit";
                    else hname = diff;
                    hname = hname + "_" + part + "_" + var + "_" + fits;
                    
                    hist_diff_part_var[fits][diff][part][var] = new TH1D(
                            hname.c_str(), hname.c_str(), 100, lbound, rbound);

                }
            }
        }
    }
}

void convert_vecs_to_car()
{
    //Given that the data vecs have already been filled by Analyzer,
    //Convert them into Cartesian for the data varcar vecs
    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = dataset.begin(); d != dataset.end(); ++d){
            const string data = *d;
            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p;
                for (unsigned int nEvent = 0; nEvent < vec_data_part_var[fits][data][part]["Pt"].size(); nEvent++){
                    //Convert to Cartesian
                    double pt = vec_data_part_var[fits][data][part]["Pt"].at(nEvent);
                    double phi = vec_data_part_var[fits][data][part]["Phi"].at(nEvent);
                    double eta = vec_data_part_var[fits][data][part]["Eta"].at(nEvent);
                    double mass = vec_data_part_var[fits][data][part]["M"].at(nEvent);
                    TLorentzVector V;
                    V.SetPtEtaPhiM(pt, eta, phi, mass);
                    double px = V.Px();
                    double py = V.Py();
                    double pz = V.Pz();
                    double energy = V.E();
 
                    vec_data_part_varcar[fits][data][part]["Px"].push_back(px);
                    vec_data_part_varcar[fits][data][part]["Py"].push_back(py);
                    vec_data_part_varcar[fits][data][part]["Pz"].push_back(pz);
                    vec_data_part_varcar[fits][data][part]["E"].push_back(energy);
                }
            }
        }
    }
}

void fill_hists()
{
    std::cout << "Fill hists" << std::endl;
    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
            const string diff = *d;

            for (auto p = particles.begin(); p != particles.end(); ++p){
                const string part = *p; 
                for (auto v = variablesH.begin(); v != variablesH.end(); ++v){
                    const string var = *v;
                    //Loop over events in vector
                    for (unsigned int nEvent = 0; nEvent < vec_diff_part_var[fits][diff][part]["Pt"].size(); nEvent++){
                        //For residual, cylindrical and M
                        if (var == "Pt" || var == "Phi" || var == "Eta" || var == "M") {
                            hist_diff_part_var[fits][diff][part][var]->Fill(vec_diff_part_var[fits][diff][part][var].at(nEvent));
                        } else if (var == "Px" or var == "Py" or var == "Pz" or var == "E") {
                            //For residual, Cartesian and E
                            double comparebase = 0;
                            double compare = 0;
                            if (diff == "smeared_gen"){
                                comparebase = (vec_data_part_varcar[fits]["gen"][part][var]).at(nEvent);
                                compare = vec_data_part_varcar[fits]["smeared"][part][var].at(nEvent);
                            } else if (diff == "best_gen"){
                                comparebase = vec_data_part_varcar[fits]["gen"][part][var].at(nEvent);
                                compare = vec_data_part_varcar[fits]["best"][part][var].at(nEvent);
                            } else if (diff == "best_smeared"){
                                comparebase = vec_data_part_varcar[fits]["smeared"][part][var].at(nEvent);
                                compare = vec_data_part_varcar[fits]["best"][part][var].at(nEvent);
                            }
                            //Fill hists
                            hist_diff_part_var[fits][diff][part][var]->Fill(compare - comparebase);
                            if (var == "Px"){
                                //std::cout << compare-comparebase << std::endl;
                            }
                        } else if (var == "Pt_" or var == "Phi_" or var == "Eta_" or var == "M_") {
                            //For resolutions, cylindrical and M
                            double comparebase = 0;
                            double compare = 0;
                            //Remove _ from end of variable name
                            const string var1 = var.substr(0, var.size()-1);

                            if (diff == "smeared_gen"){
                                comparebase = vec_data_part_var[fits]["gen"][part][var1].at(nEvent);
                                compare = vec_data_part_var[fits]["smeared"][part][var1].at(nEvent);
                            } else if (diff == "best_gen"){
                                comparebase = vec_data_part_var[fits]["gen"][part][var1].at(nEvent);
                                compare = vec_data_part_var[fits]["best"][part][var1].at(nEvent);
                            } else if (diff == "best_smeared"){
                                comparebase = vec_data_part_var[fits]["smeared"][part][var1].at(nEvent);
                                compare = vec_data_part_var[fits]["best"][part][var1].at(nEvent);
                            }
                            hist_diff_part_var[fits][diff][part][var]->Fill( (compare - comparebase)/comparebase );
                        } else if (var == "Px_" or var == "Py_" or var == "Pz_" or var == "E_") {
                            //For resolutions, Cartesian and E
                            double comparebase = 0;
                            double compare = 0;
                            //Remove _ from end of variable name
                            const string var1 = var.substr(0, var.size()-1);

                            if (diff == "smeared_gen"){
                                comparebase = vec_data_part_varcar[fits]["gen"][part][var1].at(nEvent);
                                compare = vec_data_part_varcar[fits]["smeared"][part][var1].at(nEvent);
                            } else if (diff == "best_gen"){
                                comparebase = vec_data_part_varcar[fits]["gen"][part][var1].at(nEvent);
                                compare = vec_data_part_varcar[fits]["best"][part][var1].at(nEvent);
                            } else if (diff == "best_smeared"){
                                comparebase = vec_data_part_varcar[fits]["smeared"][part][var1].at(nEvent);
                                compare = vec_data_part_varcar[fits]["best"][part][var1].at(nEvent);
                            }
                            hist_diff_part_var[fits][diff][part][var]->Fill( (compare - comparebase)/comparebase );
                        }
                    }
                }
            }
        }
    }
}

void declare_canvases()
{
    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto p = particles.begin(); p != particles.end(); ++p){
            const string part = *p;
            for (auto v = variablesH.begin(); v != variablesH.end(); ++v){
                const string var = *v;
                const string cname = "c_" + fits + "_" + part + "_" + var;
                canvs_part_var[fits][part][var] = new TCanvas(cname.c_str(), cname.c_str(), 700, 700);
            }
        }
    }
}

void plot_hists(const string dir)
{
    cout << "Drawing histograms..." << endl;
    
    //Line Colours
    int mc1 = 5;
    int mc2 = 4;
    
    int iter = 0; //To keep track of first and last plot, when combining into single PDF file
    
    for (auto fs = fitstatus.begin(); fs != fitstatus.end(); ++fs){
        const string fits = *fs;
        for (auto p = particles.begin(); p != particles.end(); ++p){
            const string part = *p;
            for (auto v = variablesH.begin(); v != variablesH.end(); ++v){
                const string var = *v;
                
                string unit = "";
                if (var == "Pt" or var == "Px" or var == "Py" or
                    var == "Pz" or var == "M" or var == "E") {
                    unit = "(GeV)";
                }
                TH1D *h_smeared = hist_diff_part_var[fits]["smeared_gen"][part][var];
                TH1D *h_best = hist_diff_part_var[fits]["best_gen"][part][var];
                TCanvas *canv = canvs_part_var[fits][part][var];

                h_smeared->SetFillColor(mc1);
                h_smeared->SetLineColor(mc1);
                h_best->SetLineColor(mc2);

                if (var == "Pt" or var == "Px" or var == "Py" or
                    var == "Phi" or var == "Pz" or var == "M" or
                    var == "E" or var == "Eta") {
                    h_smeared->GetXaxis()->SetTitle(
                        (var + " Residual " + unit).c_str());
                } else {
                    h_smeared->GetXaxis()->SetTitle(
                        (var + " Resolution " + unit).c_str());
                }
                h_smeared->GetYaxis()->SetTitle("Events");
                h_smeared->GetYaxis()->SetTitleOffset(1.4);
                //h_smeared->GetYaxis()->SetTitleSize(.1);
                //h_smeared->GetXaxis()->SetTitleFont(12);
                //h_smeared->SetTitle((fits + "_" + part + "_" + var).c_str());
                TGaxis *myY = (TGaxis*)h_smeared->GetXaxis();
                myY->SetMaxDigits(3);
                h_smeared->GetXaxis()->SetNoExponent(kTRUE);
                h_smeared->SetTitle((part + " " + var +  " (" + fits + ")").c_str());
                h_smeared->SetMaximum(
                    1.05*max(h_smeared->GetMaximum(), h_best->GetMaximum()) );
                
                

                canv->cd();

                h_smeared->Draw("HIST");
                h_best->Draw("SAMES");

                moveStatsBox(h_smeared);

                canv->Write();
                canv->Print((dir + "c_" + fits + "_" + part + "_" + var + ".pdf").c_str());

                if (iter == 0) {
                    canv->Print( (dir + "plots.pdf[").c_str() );
                    canv->Print( (dir + "plots.pdf").c_str() );
                } else if (iter ==
                           (int)(fitstatus.size()) * (int)(particles.size()) * (int)(variables.size()) - 1) {
                    canv->Print( (dir + "plots.pdf").c_str() );
                    canv->Print( (dir + "plots.pdf]").c_str() );
                } else {
                    canv->Print( (dir + "plots.pdf").c_str() );
                }
                canv->ls();

                ++iter;
            }
        }
    }

}

void read_files(const string bigpath)
{    
    initialise_vecs();

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
                    //cout << "opening file " << (bigpath + diff + "_" + part + "_" + var + "_" + fits + ".txt").c_str() << endl;
                    while (getline(infile, line)){
                        //cout << line << endl;
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
    string bigpath = "/afs/cern.ch/work/s/shtan/private/topreco_20161213/20170213/nonexisting_folder/" + path_suffix + "/";

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

                    file_diff_part_var[fits][diff][part][var] = tempfile;
                    file_diff_part_var[fits][diff][part][var]->open (
                            (bigpath + diff + "_" + part + "_" + var + "_" + fits + ".txt").c_str() );

                    *file_diff_part_var[fits][diff][part][var] << "ha" << endl;
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
    //No longer used
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

    TCanvas* canvas = new TCanvas("canvas");
    TGraph *g = new TGraph(n, arrx, arry);
    g->SetMarkerStyle(2);

    g->Draw("ap");
    g->GetXaxis()->SetTitle("Log(outer_edm)");
    g->GetYaxis()->SetTitle("Log(diffchi2_total)");
    //g->GetXaxis()->SetRangeUser(0,1);
    //g->GetYaxis()->SetRangeUser(0,20000);
    g->SetTitle("Converged");
    //g->Draw("ap");
    //canvas->Update();

    canvas->Print(filename.c_str());


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
        }
    }

    TCanvas* canvas = new TCanvas("canvas");
    TGraph *g = new TGraph(n, arrx, arry);
    g->SetMarkerStyle(2);

    g->Draw("ap");
    g->GetXaxis()->SetTitle("Log(outer_edm)");
    g->GetYaxis()->SetTitle("best_gen_diffchi2_total - smeared_gen_diffchi2_total");
    //g->GetXaxis()->SetRangeUser(0,1);
    //g->GetYaxis()->SetRangeUser(0,20000);
    g->SetTitle("Converged");
    g->Draw("ap");
    canvas->Update();

    //canvas->Print(filename.c_str());
    canvas->Print("test.pdf");

}

void plot()
{

    string plotpath = "./scatterplots/20170220/goodbad/";
    string ypart, yvar;

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

void moveStatsBox(TH1D *hist)
{
    //Move stats box in plot
    cout << "in stats box func" << endl;
    gPad->Update();
    cout << "after gpad update" << endl;
    TPaveStats *s = (TPaveStats *)hist->FindObject("stats");
    cout << "1" << endl;

    float y1 = s->GetY1NDC();
    float y2 = s->GetY2NDC();
    cout << "2" << endl;
    s->SetY1NDC(y1 - (y2 - y1));
    s->SetY2NDC(y2 - (y2 - y1));
    cout << "3" << endl;
}

#endif
