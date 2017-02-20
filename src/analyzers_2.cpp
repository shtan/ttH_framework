#ifndef analyzers_2_cpp
#define analyzers_2_cpp

#include "analyzers_2.h"
#include <iostream>
#include <cmath>
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"

#define PI acos(-1.0)

namespace ex1
{
/* 
 * analyzer2
 */
analyzer2::analyzer2()
{
    // adjust your smearer's resolution
    smr.j_res.pT = 0.1;
    smr.j_res.phi = 0.01;
    smr.j_res.eta = 0.01;

    smr.l_res.pT = 0.01;
    smr.l_res.phi = 0.001;
    smr.l_res.eta = 0.001;
    
    params.m_t = 173;
    params.m_W = 80.4;
    params.SD_m_t = 2.0;
    params.SD_m_W = 2.085;

    Setup_Maps();
}

analyzer2::~analyzer2()
{
    for (auto i : ps)
        delete i.second;
}

/*
 * analyzer2::analz() block
 * 
 */

void analyzer2::Setup_Maps()
{
/*    for (auto p = particles.begin(); p != particles.end(); ++p){
        const string part = *p;
        for (auto v = variables.begin(); v != variables.end(); ++v){
            const string var = *v;
            best_gen[part][var] = 0.0;
            smeared_gen[part][var] = 0.0;
        }
    }*/
    for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
        const string diff = *d;
        for (auto p = particles.begin(); p != particles.end(); ++p){
            const string part = *p;
            for (auto v = variables.begin(); v != variables.end(); ++v){
                const string var = *v;
                diff_part_var[diff][part][var] = 0.0;
            }
        }
    }

    for (auto d = dataset.begin(); d != dataset.end(); ++d){
        const string data = *d;
        for (auto p = particles.begin(); p != particles.end(); ++p){
            const string part = *p;
            for (auto v = variables.begin(); v != variables.end(); ++v){
                const string var = *v;
                data_part_var[data][part][var] = 0.0;
            }
        }
    }

    for (auto s = singleints.begin(); s != singleints.end(); ++s){
        const string single = *s;
        singleint[single] = -999;
    }

    for (auto s = singledoubles.begin(); s != singledoubles.end(); ++s){
        const string single = *s;
        singledouble[single] = -999;
    }

    for (auto c = chi2s.begin(); c != chi2s.end(); ++c){
        const string chi = *c;
        chisquares[chi] = 1.e99;
    }

    for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
        const string diff = *d;
        for (auto v = diffvals.begin(); v != diffvals.end(); ++v){
            const string val = *v;
            diff_diffvals[diff][val] = 0.0;
        }
    }


}

inline void CartTXYZ_to_cyl(const double *pTXYZ, double *pcyl)
{
    const double pT = sqrt(pTXYZ[1] * pTXYZ[1] + pTXYZ[2] * pTXYZ[2]);
    const double eta = asinh(pTXYZ[3] / pT);
    double phi = asin(pTXYZ[2] / pT);
    if (pTXYZ[1] < 0)
        phi = PI - phi;
    pcyl[3] = pTXYZ[0];
    pcyl[0] = pT;
    pcyl[1] = phi;
    pcyl[2] = eta;
}

inline void P4_resize(vector<double *> &p, const unsigned int sz)
{
    const auto old_sz = p.size();
    if (old_sz < sz) {
        p.reserve(sz);
        for (auto i = old_sz; i < sz; ++i)
            p.push_back(new double[4]);

        return;
    }
    if (old_sz > sz) {
        for (auto i = (old_sz - 1); i >= sz; --i)
            delete p[i];
        p.resize(sz);

        return;
    }
}

inline void Fill_input_top(const pvec &t, bool &lepton, double b[4],
                           double d1[4], double d2[4], double &met_px, double &met_py,
                           double neu[4])
{
    for (const auto id_p : t) {
        const auto &id = id_p.first;
        if (id == 5 || id == -5) {
            CartTXYZ_to_cyl(id_p.second, b);
            met_px -= id_p.second[1];
            met_py -= id_p.second[2];
            continue;
        }
        if (id == 11 || id == -11 || id == 13 || id == -13) {
            lepton = true;
            CartTXYZ_to_cyl(id_p.second, d1);
            met_px -= id_p.second[1];
            met_py -= id_p.second[2];
            continue;
        }
        /*if (id == 12 || id == -12 || id == 14 || id == -14) {
            CartTXYZ_to_cyl(id_p.second, d2);
            continue;
        }*/
        if (id == 12 || id == -12 || id == 14 || id == -14) {
            CartTXYZ_to_cyl(id_p.second, neu);
            continue;
        }
        if (id >= -4 && id <= 4 && id != 0) {
            if (d1[0] == 0) {    // check
                CartTXYZ_to_cyl(id_p.second, d1);
                met_px -= id_p.second[1];
                met_py -= id_p.second[2];
                continue;
            }
            CartTXYZ_to_cyl(id_p.second, d2);
            met_px -= id_p.second[1];
            met_py -= id_p.second[2];
            continue;
        }
    }
}

inline void Fill_input(const ttbarX &in, recoc::input_x<4> &out, recoc::input &outbig)
{
    out.t1_lep = false;
    out.t2_lep = false;
    out.d11[0] = 0;   // check
    out.d21[0] = 0;   // check
    
    out.d12[0] = 0;
    out.d12[1] = 0;
    out.d12[2] = 0;
    out.d12[3] = 0;
    out.d22[0] = 0;
    out.d22[1] = 0;
    out.d22[2] = 0;
    out.d22[3] = 0;

    double met_px = 0;
    double met_py = 0;

    Fill_input_top(in.p_t1, out.t1_lep, out.b1, out.d11, out.d12, met_px, met_py, outbig.neu_p.n1);
    Fill_input_top(in.p_t2, out.t2_lep, out.b2, out.d21, out.d22, met_px, met_py, outbig.neu_p.n2);
    
    // link over stuff
    unsigned int sz;
    if (in.p_X.size() >= 2)
        sz = in.p_X.size() - 2;
    else
        sz = 0;
    P4_resize(out.p_others, sz);
    
    auto it1 = out.p_others.begin();
    for (const auto id_p : in.p_X) {
        const auto &id = id_p.first;
        if (id == 5) {
            CartTXYZ_to_cyl(id_p.second, out.bH1);
            met_px -= id_p.second[1];
            met_py -= id_p.second[2];
            continue;
        }
        if (id == -5) {
            CartTXYZ_to_cyl(id_p.second, out.bH2);
            met_px -= id_p.second[1];
            met_py -= id_p.second[2];
            continue;
        }
        
        // "exotics" go in here
        CartTXYZ_to_cyl(id_p.second, *it1);
        met_px -= id_p.second[1];
        met_py -= id_p.second[2];
        ++it1;
    }

    outbig.MET_px = met_px;
    outbig.MET_py = met_py;

}

inline void Fill_SDs_(const double p[4], const double unc[3], double SD[3])
{
    SD[0] = unc[0] * p[0];
    SD[1] = unc[1] * p[1];
    SD[2] = unc[2] * p[2];
}

void Fill_SDs(recoc::input &in)
{
    cout << "inside Fill_SDs" << endl;
    auto &p = in.p;
    auto &SD = in.SD;
    const double j_res[3] = {0.1, 0.01, 0.01};  // jet: pT, phi, eta
    const double l_res[3] = {0.01, 0.001, 0.001};   // lepton: pT, phi, eta
    
    Fill_SDs_(p.b1, j_res, SD.b1);
    Fill_SDs_(p.b2, j_res, SD.b2);
    Fill_SDs_(p.bH1, j_res, SD.bH1);
    Fill_SDs_(p.bH2, j_res, SD.bH2);
    if (p.t1_lep)
        Fill_SDs_(p.d11, l_res, SD.d11);
    else {
        Fill_SDs_(p.d11, j_res, SD.d11);
        Fill_SDs_(p.d12, j_res, SD.d12);
    }
    if (p.t2_lep)
        Fill_SDs_(p.d21, l_res, SD.d21);
    else {
        Fill_SDs_(p.d21, j_res, SD.d21);
        Fill_SDs_(p.d22, j_res, SD.d22);
    }
    
//    if (p.p_others.size() == 0)
//        return;

    auto it1 = p.p_others.begin();
    auto it2 = SD.p_others.begin();
    const auto end = p.p_others.end();
    while (it1 != end) {
        Fill_SDs_(*it1, j_res, *it2);
    }

    cout << "after original" << endl;
    double metSD2[3] = {0., 0., 0.};

    for (int i = 0; i<3; i++){
        for (auto it = SD.p_others.begin(); it < SD.p_others.end(); ++it ){
            metSD2[i] += pow( (*it)[i], 2);
        }

        metSD2[i] += pow( SD.b1[i], 2);
        cout << "SD.b1 = " << SD.b1[0] << " " << SD.b1[1] << " " << SD.b1[2] << endl; 
        metSD2[i] += pow( SD.b2[i], 2);
        metSD2[i] += pow( SD.bH1[i], 2);
        metSD2[i] += pow( SD.bH2[i], 2);
        metSD2[i] += pow( SD.d11[i], 2);
        if (!p.t1_lep){
            metSD2[i] += pow( SD.d12[i], 2);
        }
        metSD2[i] += pow( SD.d21[i], 2);
        if (!p.t2_lep){
            metSD2[i] += pow( SD.d22[i], 2);
        }
    }

    double metSD[3] = { pow(metSD2[0], 0.5), pow(metSD2[1], 0.5), pow(metSD2[2], 0.5) };
    if (p.t1_lep){
        memcpy(SD.d12, metSD, sizeof(SD.d12));
    }
    if (p.t2_lep){
        memcpy(SD.d22, metSD, sizeof(SD.d22));
    }

    cout << "SD.d12 = " << SD.d12[0] << " " << SD.d12[1] << " " << SD.d12[2] << endl;

    //return;
}

void analyzer2::smear_only(const pvec &p)
{
    smr.smear(p, ps);
}

void analyzer2::analz(const pvec &p, const movec &moth_ID, long unsigned int &event_num,
                        vdmap4 &vec_diff_part_var, vdmap4 &vec_data_part_var, vimap2 &vec_singleint,
                        vdmap2 &vec_singledouble, vdmap2 &vec_chisquares, vdmap3 &vec_diff_diffvals)
/*void analyzer2::analz(const pvec &p, const movec &moth_ID, long unsigned int &event_num, fmap2 &outfiles_best_gen_all, fmap2 &outfiles_smeared_gen_all,
                        fmap2 &outfiles_best_gen_converged, fmap2 &outfiles_smeared_gen_converged,
                        fmap2 &outfiles_best_gen_failed, fmap2 &outfiles_smeared_gen_failed,
                        fmap2 &outfiles_best_all, fmap2 &outfiles_smeared_all, fmap2 &outfiles_gen_all,
                        ofstream &outfile_inner_status_all, ofstream &outfile_outer_status_all, ofstream &outfile_event_number_all,
                        ofstream &outfile_inner_edm_all, ofstream &outfile_outer_edm_all,
                        fmap2 &outfiles_best_converged, fmap2 &outfiles_smeared_converged, fmap2 &outfiles_gen_converged,
                        ofstream &outfile_inner_status_converged, ofstream &outfile_outer_status_converged, ofstream &outfile_event_number_converged,
                        ofstream &outfile_inner_edm_converged, ofstream &outfile_outer_edm_converged,
                        fmap2 &outfiles_best_failed, fmap2 &outfiles_smeared_failed, fmap2 &outfiles_gen_failed,
                        ofstream &outfile_inner_status_failed, ofstream &outfile_outer_status_failed, ofstream &outfile_event_number_failed,
                        ofstream &outfile_inner_edm_failed, ofstream &outfile_outer_edm_failed,
                        ofstream &outfile_smeared_gen_diffchi2_all, ofstream &outfile_best_gen_diffchi2_all, ofstream &outfile_best_smeared_diffchi2_all,
                        ofstream &outfile_smeared_gen_diffchi2_converged, ofstream &outfile_best_gen_diffchi2_converged, ofstream &outfile_best_smeared_diffchi2_converged,
                        ofstream &outfile_smeared_gen_diffchi2_failed, ofstream &outfile_best_gen_diffchi2_failed, ofstream &outfile_best_smeared_diffchi2_failed
                        )*/
{
    smr.smear(p, ps);
    ttbarX ttX = ident.identify(ex1::ttH_SL_bx, p, moth_ID);
    ttbarX ttX_smr = ident.identify(ex1::ttH_SL_bx, ps, moth_ID);
    
    cout << endl;
    cout << "ttX" << endl;
    ttX.Print_contents();
     cout << endl;
    cout << "ttX_smr" << endl;
    ttX_smr.Print_contents();
    cout << endl;
    
    Fill_input(ttX, generated.p, generated);
    Fill_input(ttX_smr, in_2_RC.p, in_2_RC);
    //Fill_input(ttX, in_2_RC.p, in_2_RC);

    Fill_SDs(generated);
    Fill_SDs(in_2_RC);
    cout<<endl;
cout <<"generated:" << endl;
    recoc::Print(generated);
    cout<<endl;
    cout << "smeared:" << endl;
    recoc::Print(in_2_RC);
    cout<<endl;
    //recoc::output result = reco_C.reco(generated, params);
    recoc::output result = reco_C.reco(in_2_RC, params);
cout <<"blahh in_2_RC" << endl;
    //recoc::Print(in_2_RC);
cout << "result" << endl;
     recoc::Print(result);
    //recoc::Print_diff(result, in_2_RC);
    cout << "DIFFERENCE BETWEEN SMEARED AND GEN" << endl;
    recoc::Print_diff(in_2_RC, generated);
    cout << "DIFFERENCE BETWEEN FITTED AND GEN" << endl;
    recoc::Print_diff(result, generated);
   
    //outfile << result.inner_min_status << " " << result.outer_min_status << endl;

    recoc::output genforcompare;
    recoc::output in_2_RC_forcompare;
    //reco_C.input_to_output (generated, genforcompare);
    //reco_C.input_to_output (in_2_RC, in_2_RC_forcompare);

    genforcompare.p = generated.p;
    in_2_RC_forcompare.p = in_2_RC.p;

    for (int i = 0; i<4; ++i){
        genforcompare.p.d12[i] = generated.neu_p.n1[i];
        cout << "generated neutrinos:" << endl;
        cout << genforcompare.p.d12[i] << endl;
    }
    reco_C.met_to_neutrino( in_2_RC.MET_px, in_2_RC.MET_py, in_2_RC_forcompare.p.d12 );
    cout << "met_px = " << in_2_RC.MET_px << endl;
    cout << "met_py = " << in_2_RC.MET_py << endl;
    for (int i = 0; i<4; ++i){
        cout << in_2_RC_forcompare.p.d12[i] << endl;
    }

    reco_C.daughter_to_parents(genforcompare.p, genforcompare.parents_p);
    reco_C.daughter_to_parents(in_2_RC_forcompare.p, in_2_RC_forcompare.parents_p);

    genforcompare.SD = generated.SD;
    in_2_RC_forcompare.SD = in_2_RC.SD;

    cout << "result - genforcompare" << endl;
    //calc_diff(result, genforcompare, best_gen);
    calc_diff(result, genforcompare, diff_part_var["best_gen"]);
    cout << "in2rc - genforcompare" << endl;
    //calc_diff(in_2_RC_forcompare, genforcompare, smeared_gen);
    //calc_diff(result, in_2_RC_forcompare, best_smeared);
    calc_diff(in_2_RC_forcompare, genforcompare, diff_part_var["smeared_gen"]);
    calc_diff(result, in_2_RC_forcompare, diff_part_var["best_smeared"]);

    dmap2_converter(result, data_part_var["best"]);
    dmap2_converter(in_2_RC_forcompare, data_part_var["smeared"]);
    dmap2_converter(genforcompare, data_part_var["gen"]);

    diff_diffvals["smeared_gen"]["diff_chi2_total"] = diff_chi2(diff_part_var["smeared_gen"], data_part_var["gen"], genforcompare, "total");
    diff_diffvals["best_gen"]["diff_chi2_total"] = diff_chi2(diff_part_var["best_gen"], data_part_var["gen"], genforcompare, "total");
    diff_diffvals["best_smeared"]["diff_chi2_total"] = diff_chi2(diff_part_var["best_smeared"], data_part_var["smeared"], in_2_RC_forcompare, "total");

    diff_diffvals["smeared_gen"]["diff_chi2_measurable"] = diff_chi2(diff_part_var["smeared_gen"], data_part_var["gen"], genforcompare, "measurable");
    diff_diffvals["best_gen"]["diff_chi2_measurable"] = diff_chi2(diff_part_var["best_gen"], data_part_var["gen"], genforcompare, "measurable");
    diff_diffvals["best_smeared"]["diff_chi2_measurable"] = diff_chi2(diff_part_var["best_smeared"], data_part_var["smeared"], in_2_RC_forcompare, "measurable");

    diff_diffvals["smeared_gen"]["diff_chi2_masses"] = diff_chi2(diff_part_var["smeared_gen"], data_part_var["gen"], genforcompare, "masses");
    diff_diffvals["best_gen"]["diff_chi2_masses"] = diff_chi2(diff_part_var["best_gen"], data_part_var["gen"], genforcompare, "masses");
    diff_diffvals["best_smeared"]["diff_chi2_masses"] = diff_chi2(diff_part_var["best_smeared"], data_part_var["smeared"], in_2_RC_forcompare, "masses");

    diff_diffvals["smeared_gen"]["diff_chi2_neutrino"] = diff_chi2(diff_part_var["smeared_gen"], data_part_var["gen"], genforcompare, "neutrino");
    diff_diffvals["best_gen"]["diff_chi2_neutrino"] = diff_chi2(diff_part_var["best_gen"], data_part_var["gen"], genforcompare, "neutrino");
    diff_diffvals["best_smeared"]["diff_chi2_neutrino"] = diff_chi2(diff_part_var["best_smeared"], data_part_var["smeared"], in_2_RC_forcompare, "neutrino");

    diff_diffvals["smeared_gen"]["diff_chi2_measurableplusmasses"] = diff_chi2(diff_part_var["smeared_gen"], data_part_var["gen"], genforcompare, "measurableplusmasses");
    diff_diffvals["best_gen"]["diff_chi2_measurableplusmasses"] = diff_chi2(diff_part_var["best_gen"], data_part_var["gen"], genforcompare, "measurableplusmasses");
    diff_diffvals["best_smeared"]["diff_chi2_measurableplusmasses"] = diff_chi2(diff_part_var["best_smeared"], data_part_var["smeared"], in_2_RC_forcompare, "measurableplusmasses");

    get_chi2s(result, chisquares);

    push_big("all", result, event_num, vec_diff_part_var, vec_data_part_var, vec_singleint,
                vec_singledouble, vec_chisquares, vec_diff_diffvals);

    if ( (result.inner_min_status == 0 or result.inner_min_status == 1)
      and (result.outer_min_status == 0 or result.outer_min_status == 1) ){
        push_big("converged", result, event_num, vec_diff_part_var, vec_data_part_var, vec_singleint,
                vec_singledouble, vec_chisquares, vec_diff_diffvals);
    } else {
        push_big("failed", result, event_num, vec_diff_part_var, vec_data_part_var, vec_singleint,
                vec_singledouble, vec_chisquares, vec_diff_diffvals);
    }

/*    write_big("all", result, event_num, file_diff_part_var, file_data_part_var, file_singleint,
                file_singledouble, file_chisquares, file_diff_diffvals);

    if ( (result.inner_min_status == 0 or result.inner_min_status == 1)
      and (result.outer_min_status == 0 or result.outer_min_status == 1) ){
        write_big("converged", result, event_num, file_diff_part_var, file_data_part_var, file_singleint,
                file_singledouble, file_chisquares, file_diff_diffvals);
    } else {
        write_big("failed", result, event_num, file_diff_part_var, file_data_part_var, file_singleint,
                file_singledouble, file_chisquares, file_diff_diffvals);
    }*/


/*    cout << "out here!" << endl;
    write_diff(best_gen, outfiles_best_gen_all);
    write_diff(smeared_gen, outfiles_smeared_gen_all);

    write_diff(best, outfiles_best_all);
    write_diff(smeared, outfiles_smeared_all);
    write_diff(gen, outfiles_gen_all);
    write_int(result.inner_min_status, outfile_inner_status_all);
    write_int(result.outer_min_status, outfile_outer_status_all);
    write_int(result.inner_edm, outfile_inner_edm_all);
    write_int(result.outer_edm, outfile_outer_edm_all);
    write_int(event_num, outfile_event_number_all);
    cout << "before write new" << endl;
    write_int(smeared_gen_diffchi2, outfile_smeared_gen_diffchi2_all);
    cout << "after first write new" << endl;
    write_int(best_gen_diffchi2, outfile_best_gen_diffchi2_all);
    write_int(best_smeared_diffchi2, outfile_best_smeared_diffchi2_all);
    cout << "after third write new" << endl;

    if ( (result.inner_min_status == 0 or result.inner_min_status == 1)
      and (result.outer_min_status == 0 or result.outer_min_status == 1) ){
        write_diff(best_gen, outfiles_best_gen_converged);
        write_diff(smeared_gen, outfiles_smeared_gen_converged);
        write_diff(best, outfiles_best_converged);
        write_diff(smeared, outfiles_smeared_converged);
        write_diff(gen, outfiles_gen_converged);
        write_int(result.inner_min_status, outfile_inner_status_converged);
        write_int(result.outer_min_status, outfile_outer_status_converged);
        write_int(result.inner_edm, outfile_inner_edm_converged);
        write_int(result.outer_edm, outfile_outer_edm_converged);
        write_int(event_num, outfile_event_number_converged);
        write_int(smeared_gen_diffchi2, outfile_smeared_gen_diffchi2_converged);
        write_int(best_gen_diffchi2, outfile_best_gen_diffchi2_converged);
        write_int(best_smeared_diffchi2, outfile_best_smeared_diffchi2_converged);
    } else {
        write_diff(best_gen, outfiles_best_gen_failed);
        write_diff(smeared_gen, outfiles_smeared_gen_failed);
        write_diff(best, outfiles_best_failed);
        write_diff(smeared, outfiles_smeared_failed);
        write_diff(gen, outfiles_gen_failed);
        write_int(result.inner_min_status, outfile_inner_status_failed);
        write_int(result.outer_min_status, outfile_outer_status_failed);
        write_int(result.inner_edm, outfile_inner_edm_failed);
        write_int(result.outer_edm, outfile_outer_edm_failed);
        write_int(event_num, outfile_event_number_failed);
        write_int(smeared_gen_diffchi2, outfile_smeared_gen_diffchi2_failed);
        write_int(best_gen_diffchi2, outfile_best_gen_diffchi2_failed);
        write_int(best_smeared_diffchi2, outfile_best_smeared_diffchi2_failed);
    }

   cout << "after all write" << endl; */


    for (auto array : result.p.p_others)
        delete array;
    cout << "after delete" << endl;
}

/*void analyzer2::input_to_output(const recoc::input &in, recoc::output &out)
{
    out.p = in.p;
    recoc::daughter_to_parents(in.p, out.parents_p);
}*/



void analyzer2::calc_diff(const recoc::output &out1, const recoc::output &out2, dmap2 &diff)
{
    //Top 1
    one_diff( out1.p.b1, out2.p.b1, diff, "Bottom_1" );
    one_diff( out1.p.d11, out2.p.d11, diff, "Wd11" );
    one_diff( out1.p.d12, out2.p.d12, diff, "Wd12" );
    one_diff( out1.parents_p.w1, out2.parents_p.w1, diff, "W1" );
    one_diff( out1.parents_p.t1, out2.parents_p.t1, diff, "Top_1" );
    //Top 2
    one_diff( out1.p.b2, out2.p.b2, diff, "Bottom_2" );
    one_diff( out1.p.d21, out2.p.d21, diff, "Wd21" );
    one_diff( out1.p.d22, out2.p.d22, diff, "Wd22" );
    one_diff( out1.parents_p.w2, out2.parents_p.w2, diff, "W2" );
    one_diff( out1.parents_p.t2, out2.parents_p.t2, diff, "Top_2" );
    //Higgs
    one_diff( out1.p.bH1, out2.p.bH1, diff, "b1_from_H" );
    one_diff( out1.p.bH2, out2.p.bH2, diff, "b2_from_H" );
    one_diff( out1.parents_p.h, out2.parents_p.h, diff, "Higgs" );
}

void analyzer2::dmap2_converter(const recoc::output &out, dmap2 &dmapp)
{
    //Top1
    single_converter( out.p.b1, dmapp, "Bottom_1" );
    single_converter( out.p.d11, dmapp, "Wd11" );
    single_converter( out.p.d12, dmapp, "Wd12" );
    single_converter( out.parents_p.w1, dmapp, "W1" );
    single_converter( out.parents_p.t1, dmapp, "Top_1" );
    //Top 2
    single_converter( out.p.b2, dmapp, "Bottom_2" );
    single_converter( out.p.d21, dmapp, "Wd21" );
    single_converter( out.p.d22, dmapp, "Wd22" );
    single_converter( out.parents_p.w2, dmapp, "W2" );
    single_converter( out.parents_p.t2, dmapp, "Top_2" );
    //Higgs
    single_converter( out.p.bH1, dmapp, "b1_from_H" );
    single_converter( out.p.bH2, dmapp, "b2_from_H" );
    single_converter( out.parents_p.h, dmapp, "Higgs" );
}


void analyzer2::principal_angle(double &theta)
{
    cout << "theta = " << theta << endl;
    double pii = 3.14159265359;
    while (theta > pii){
        theta -= 2*pii;
    }
    while (theta < -pii){
        theta += 2*pii;
    }
    cout << "new theta = " << theta << endl;
}

void analyzer2::single_converter(const double vec[4], dmap2 &dmapp, string partname)
{
    dmapp[partname]["Pt"] = vec[0];

    double phi = vec[1];
    double eta = vec[2];
    principal_angle(phi);
    //principal_angle(eta);

    dmapp[partname]["Phi"] = phi;
    dmapp[partname]["Eta"] = eta;

    TLorentzVector lv;
    lv.SetPtEtaPhiE(vec[0], vec[2], vec[1], vec[3]);

    dmapp[partname]["M"] = lv.M();
    cout << "blub " << partname << " " << lv.M() << endl;

}

void analyzer2::get_chi2s(recoc::output& result, dmap1 &dmap)
{
    dmap["Bottom_1"] = result.chi2s.b1;
    dmap["Wd11"] = result.chi2s.d11;
    dmap["Wd12"] = result.chi2s.d12;
    dmap["mW1"] = result.chi2s.mw1;
    dmap["mTop_1"] = result.chi2s.mt1;
    dmap["Bottom_2"] = result.chi2s.b2;
    dmap["Wd21"] = result.chi2s.d21;
    dmap["Wd22"] = result.chi2s.d22;
    dmap["mW2"] = result.chi2s.mw2;
    dmap["mTop_2"] = result.chi2s.mt2;
    dmap["nontops"] = result.chi2s.nontop_objs;
}

double analyzer2::one_p_chi2(dmap2 &diff, string partname, const double SD[3])
{
    cout << "starting one_p_chi2" << endl;
    cout << partname << endl;
    double pt = diff[partname]["Pt"];
    double phi = diff[partname]["Phi"];
    double eta = diff[partname]["Eta"];
    double SD_pt = SD[0];
    double SD_phi = SD[1];
    double SD_eta = SD[2];
    cout << "SD" << endl;
    cout << SD_pt << endl;
    cout << SD_phi << endl;
    cout << SD_eta << endl;
    cout << "P" << endl;
    cout << pt << endl;
    cout << phi << endl;
    cout << eta << endl;

    double chi2 = pow( pt/SD_pt, 2) + pow( phi/SD_phi, 2) + pow( eta/SD_eta, 2);
    cout << "chi2 = " << chi2 << endl;

    return chi2;
}

//Define breitWigner function to calculate mTop and mW chi2
double analyzer2::breitWignerErr(const double &mass, const double &width, const double &deltaMass)
{
    cout << "in breitwigner" << endl;
    double scaledDeltaMass = deltaMass * width;
    double scaledDeltaMass2 = scaledDeltaMass * scaledDeltaMass;
    double Zscore = ROOT::Math::normal_quantile(
                0.31830988618379067154 *atan2(scaledDeltaMass2 + 2. * scaledDeltaMass * mass,
                mass * width) + 0.5, 1.0);
    return Zscore * Zscore;
}


double analyzer2::one_m_chi2(dmap2 &diff, dmap2 &control, string partname, const double width)
{
    cout << "starting one_m_chi2"<< endl;
    cout << partname << endl;
    double control_mass = control[partname]["M"];
    double mass_diff = diff[partname]["M"];
    double chi2 = breitWignerErr(control_mass, width, mass_diff);
    cout << "chi2 = " << chi2 << endl;
    return chi2;
}

double analyzer2::diff_chi2(dmap2 &diff, dmap2 &control, recoc::output &compare, string which)
{
    cout << "starting diff_chi2" << endl;
    double b1_chi2 = one_p_chi2(diff, "Bottom_1", compare.SD.b1);
    double wd11_chi2 = one_p_chi2(diff, "Wd11", compare.SD.d11);
    double wd12_chi2 = one_p_chi2(diff, "Wd12", compare.SD.d12);
    double b2_chi2 = one_p_chi2(diff, "Bottom_2", compare.SD.b2);
    double wd21_chi2 = one_p_chi2(diff, "Wd21", compare.SD.d21);
    double wd22_chi2 = one_p_chi2(diff, "Wd22", compare.SD.d22);
    double bH1_chi2 = one_p_chi2(diff, "b1_from_H", compare.SD.bH1);
    double bH2_chi2 = one_p_chi2(diff, "b2_from_H", compare.SD.bH2);
    double mW1_chi2 = one_m_chi2(diff, control, "W1", params.SD_m_W);
    double mW2_chi2 = one_m_chi2(diff, control, "W2", params.SD_m_W);
    double mTop1_chi2 = one_m_chi2(diff, control, "Top_1", params.SD_m_t);
    double mTop2_chi2 = one_m_chi2(diff, control, "Top_2", params.SD_m_t);

    double chi2 = 1.e99;
    
    if (strcmp(which.c_str(), "total") == 0){
        chi2 = b1_chi2 + wd11_chi2 + wd12_chi2 
                    + b2_chi2 + wd21_chi2 + wd22_chi2
                    + bH1_chi2 + bH2_chi2
                    + mW1_chi2 + mW2_chi2 + mTop1_chi2 + mTop2_chi2;
    } else if (strcmp(which.c_str(), "measurable") == 0){
        //version without masses and without real neutrino
        chi2 = b1_chi2 + wd11_chi2 
                    + b2_chi2 + wd21_chi2 + wd22_chi2
                    + bH1_chi2 + bH2_chi2;
    } else if (strcmp(which.c_str(), "masses") == 0){
        chi2 = mW1_chi2 + mW2_chi2 + mTop1_chi2 + mTop2_chi2;
    } else if (strcmp(which.c_str(), "neutrino") == 0){
        chi2 = wd12_chi2;
    } else if (strcmp(which.c_str(), "measurableplusmasses") == 0){
        //version without real neutrino
        chi2 = b1_chi2 + wd11_chi2 
                    + b2_chi2 + wd21_chi2 + wd22_chi2
                    + bH1_chi2 + bH2_chi2
                    + mW1_chi2 + mW2_chi2 + mTop1_chi2 + mTop2_chi2;
    }

    return chi2;
}

void analyzer2::one_diff(const double vec1[4], const double vec2[4], dmap2 &diff, string partname)
{
    diff[partname]["Pt"] = vec1[0] - vec2[0];

    double phi = vec1[1] - vec2[1];
    double eta = vec1[2] - vec2[2];
    principal_angle(phi);
    //principal_angle(eta);

    diff[partname]["Phi"] = phi;
    diff[partname]["Eta"] = eta;
    //diff[partname]["Phi"] = vec1[1] - vec2[1];
    //diff[partname]["Eta"] = vec1[2] - vec2[2];
    
    //Get masses
    TLorentzVector lv1;
    lv1.SetPtEtaPhiE(vec1[0], vec1[2], vec1[1], vec1[3]);
    TLorentzVector lv2;
    lv2.SetPtEtaPhiE(vec2[0], vec2[2], vec2[1], vec2[3]);

    diff[partname]["M"] = lv1.M() - lv2.M();

    cout << "partname = " << partname << endl;
    cout << "Pt " << vec1[0] << " " << vec2[0] << endl;
    cout << "Phi " << vec1[1] << " " << vec2[1] << endl;
    cout << "Eta " << vec1[2] << " " << vec2[2] << endl;
    cout << "Px " << lv1.Px() << " " << lv2.Px() << endl;
    cout << "Py " << lv1.Py() << " " << lv2.Py() << endl;
    cout << "Pz " << lv1.Pz() << " " << lv2.Pz() << endl;
    cout << "E " << lv1.E() << " " << lv2.E() << endl;
    cout << "lv1 mass = " << lv1.M() << endl;
    cout << "lv2 mass = " << lv2.M() << endl;

}

void analyzer2::write_big(string fit_status, recoc::output &result, long unsigned int &event_num,
                            fmap4 &file_diff_part_var, fmap4 &file_data_part_var, fmap2 &file_singleint,
                            fmap2 &file_singledouble, fmap2 &file_chisquares, fmap3 &file_diff_diffvals)
{

    for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
        const string diff = *d;
        write_diff(diff_part_var[diff], file_diff_part_var[fit_status][diff]);
    }

/*    for (auto d = dataset.begin(); d != dataset.end(); ++d){
        const string data = *d;
        write_diff(data_part_var[data], file_data_part_var[fit_status][data]);
    }*/
/*    write_diff(diff_part_var["best_gen"], file_diff_part_var[fit_status]["best_gen"];
    write_diff(diff_part_var["smeared_gen"], file_diff_part_var[fit_status]["smeared_gen"]);
    
    write_diff(data_part_var["best"], file_data_part_var[fit_status]["best"]);
    write_diff(data_part_var["smeared"], file_data_part_var[fit_status]["smeared"]);
    write_diff(data_part_var["gen"], file_data_part_var[fit_status]["gen"]);*/

    write_int(result.inner_min_status, file_singleint[fit_status]["inner_status"]);
    write_int(result.outer_min_status, file_singleint[fit_status]["outer_status"]);
    write_int(result.inner_edm, file_singledouble[fit_status]["inner_edm"]);
    write_int(result.outer_edm, file_singledouble[fit_status]["outer_edm"]);
    write_int(event_num, file_singleint[fit_status]["event_num"]);

    for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
        const string diff = *d;
        for (auto v = diffvals.begin(); v != diffvals.end(); ++v){
            const string val = *v;
            write_int(diff_diffvals[diff][val], file_diff_diffvals[fit_status][diff][val]);
        }
    }

/*    for (auto c = chi2s.begin(); c != chi2s.end(); ++c){
        const string chi = *c;
        write_int(chisquares[fits][chi]*/

/*    write_int(diff_diffvals["smeared_gen"]_diffchi2, outfile_smeared_gen_diffchi2_all);
    write_int(best_gen_diffchi2, outfile_best_gen_diffchi2_all);
    write_int(best_smeared_diffchi2, outfile_best_smeared_diffchi2_all);*/
}

void analyzer2::push_big(string fit_status, recoc::output &result, long unsigned int &event_num,
                            vdmap4 &vec_diff_part_var, vdmap4 &vec_data_part_var, vimap2 &vec_singleint,
                            vdmap2 &vec_singledouble, vdmap2 &vec_chisquares, vdmap3 &vec_diff_diffvals)
{

    for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
        const string diff = *d;
        push_diff(diff_part_var[diff], vec_diff_part_var[fit_status][diff]);
    }

    for (auto d = dataset.begin(); d != dataset.end(); ++d){
        const string data = *d;
        push_diff(data_part_var[data], vec_data_part_var[fit_status][data]);
    }
/*    write_diff(diff_part_var["best_gen"], file_diff_part_var[fit_status]["best_gen"];
    write_diff(diff_part_var["smeared_gen"], file_diff_part_var[fit_status]["smeared_gen"]);
    
    write_diff(data_part_var["best"], file_data_part_var[fit_status]["best"]);
    write_diff(data_part_var["smeared"], file_data_part_var[fit_status]["smeared"]);
    write_diff(data_part_var["gen"], file_data_part_var[fit_status]["gen"]);*/

    push_int(result.inner_min_status, vec_singleint[fit_status]["inner_status"]);
    push_int(result.outer_min_status, vec_singleint[fit_status]["outer_status"]);
    push_int(result.inner_edm, vec_singledouble[fit_status]["inner_edm"]);
    push_int(result.outer_edm, vec_singledouble[fit_status]["outer_edm"]);
    push_int(event_num, vec_singleint[fit_status]["event_number"]);

    for (auto d = datasetdiff.begin(); d != datasetdiff.end(); ++d){
        const string diff = *d;
        for (auto v = diffvals.begin(); v != diffvals.end(); ++v){
            const string val = *v;
            push_int(diff_diffvals[diff][val], vec_diff_diffvals[fit_status][diff][val]);
        }
    }

    for (auto c = chi2s.begin(); c != chi2s.end(); ++c){
        const string chi = *c;
        push_int(chisquares[chi], vec_chisquares[fit_status][chi]);
    }

/*    write_int(diff_diffvals["smeared_gen"]_diffchi2, outfile_smeared_gen_diffchi2_all);
    write_int(best_gen_diffchi2, outfile_best_gen_diffchi2_all);
    write_int(best_smeared_diffchi2, outfile_best_smeared_diffchi2_all);*/
}

void analyzer2::push_diff(dmap2 &diff, vdmap2 &vecmap)
{
    for (auto p = particles.begin(); p != particles.end(); ++p){
        const string part = *p;
        for (auto v = variables.begin(); v != variables.end(); ++v){
            const string var = *v;
            (vecmap[part][var]).push_back(diff[part][var]);
            //cout << "pushing " << part << " " << endl;
        }
    }

}

void analyzer2::push_int(int &num, vector<int> &vec)
{
    vec.push_back(num);
}

void analyzer2::push_int(long unsigned int &num, vector<int> &vec)
{
    //cout << "inside long unsigned push" << endl;
    vec.push_back(int(num));
    //cout << "unsigned vec size = " << vec.size() << endl;
}

void analyzer2::push_int(double &num, vector<double> &vec)
{
    vec.push_back(num);
}



void analyzer2::write_diff(dmap2 &diff, fmap2 &outfiles)
{
    for (auto p = particles.begin(); p != particles.end(); ++p){
        const string part = *p;
        for (auto v = variables.begin(); v != variables.end(); ++v){
            const string var = *v;
            *(outfiles[part][var]) << diff[part][var] << endl;
        }
    }

}

void analyzer2::write_int(int &num, ofstream *outfile)
{
    *outfile << num << endl;
}

void analyzer2::write_int(long unsigned int &num, ofstream *outfile)
{
    *outfile << num << endl;
}

void analyzer2::write_int(double &num, ofstream *outfile)
{
    *outfile << num << endl;
}

}
#endif
