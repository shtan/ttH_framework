#ifndef scatter_plotter_cpp
#define scatter_plotter_cpp

#include "scatter_plotter.h"
#include <iostream>
#include <cmath>

#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPDF.h>
#include <TLine.h>

namespace plotter
{
/* 
 * analyzer2
 */
plotter::plotter()
{
}

plotter::~plotter()
{
}

void plotter::good_bad(vector<double> &vecx, vector<double> &best_gen, vector<double> &smeared_gen, string filename, string title, string xaxis_label)
{
    cout << "scattering plot" << endl;
    int n = vecx.size();
    double arrx[n];
    double arry[n];
    for (int i=0; i<n; i++){
        //arrx[i] = log10(vecx.at(i)+1e-20);
        arrx[i] = vecx.at(i);
        double dif = abs(best_gen.at(i)) - abs(smeared_gen.at(i));
        arry[i] = dif;
    }

    cout << "created arrs" << endl;
    TCanvas* canvas = new TCanvas("canvas");
    TGraph *g = new TGraph(n, arrx, arry);
    g->SetMarkerStyle(2);

    cout << "getting axes" << endl;
    g->Draw("ap");
    g->GetXaxis()->SetTitle(xaxis_label.c_str());
    g->GetYaxis()->SetTitle("abs(best_gen) - abs(smeared_gen)");
    //g->GetXaxis()->SetRangeUser(0,1);
    //g->GetYaxis()->SetRangeUser(0,20000);
    g->SetTitle((title).c_str());

    canvas->Update();

    TLine *l = new TLine(canvas->GetUxmin(), 0.0, canvas->GetUxmax(), 0.0);
    l->SetLineColor(kBlue);
    l->Draw();

    cout << "before draw" << endl;
    //g->Draw("ap");
    canvas->Update();

    cout << "after draw" << endl;
    canvas->Print(filename.c_str());
    cout << "after print" << endl;

}

}
#endif
