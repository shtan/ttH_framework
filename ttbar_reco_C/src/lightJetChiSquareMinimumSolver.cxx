#ifndef lightJetChiSquareMinimumSolver_cxx
#define lightJetChiSquareMinimumSolver_cxx

#include "lightJetChiSquareMinimumSolver.h"
#include "TDecompSVD.h"
#include "TDecompLU.h"

using namespace commonstruct;

lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(
    big_struct &bigstructure, bool do3D, int &debug, nontop_system &testnontop)
    : bigstruct(bigstructure), nontops( *(bigstruct.nontops_ptr) ),
      do3D_(do3D), debug_verbosity(debug), testnontops(testnontop)
{
    if (debug_verbosity >= 2)
        cout << "Starting lightJetChiSquareMinimumSolver constructor" << endl;

    Init_data(data_);

    inverter3D_ = new TDecompLU(3);

    checkSize();
    Eval_covariance(data_);
    Eval_cov_sum(data_);

}


inline void lightJetChiSquareMinimumSolver::Init_data(recoil_minimizer_data &da)
{
    if (debug_verbosity >= 2)
        cout << "Starting init_data, n_nontop_objs = " << nontops.input.n_objs << endl;
    
    //da.n_ps = sz;
    da.dxCheck = 0.;
    da.dyCheck = 0.;
    da.dzCheck = 0.;

    da.cov_rad.reserve(nontops.input.n_objs);
    da.cov.reserve(nontops.input.n_objs);
    da.inv_sum_x_cov.reserve(nontops.input.n_objs);
    for (int i = 0; i < nontops.input.n_objs; ++i) {
        da.cov_rad.push_back(TMatrixD(3, 3));
        da.cov.push_back(TMatrixD(3, 3));
        da.inv_sum_x_cov.push_back(TMatrixD(3, 3));
    }
}


lightJetChiSquareMinimumSolver::~lightJetChiSquareMinimumSolver()
{
    delete inverter3D_;
}


void lightJetChiSquareMinimumSolver::setupEquations()
{
    checkSize();
    Eval_covariance(data_);
    Eval_cov_sum(data_);
}

void lightJetChiSquareMinimumSolver::checkSize()
{
    if (debug_verbosity >= 2)
                cout << "starting checkSize" << endl;

    if (nontops.input.jet_pt.size() != nontops.input.jet_eta.size()) {
        cout << "Unequal number of jet pts and jet pT widths!" << endl;
        cout << "there are " << nontops.input.jet_pt.size() << " jet pts and "
             << nontops.input.jet_eta.size() << " jet etas" << endl;
        return;
    }
    if (nontops.input.jet_pt.size() != nontops.input.jet_phi.size()) {
        cout << "Unequal number of jet pts and jet pT widths!" << endl;
        cout << "there are " << nontops.input.jet_pt.size() << " jet pts and "
             << nontops.input.jet_phi.size() << " jet phis" << endl;
        return;
    }
    if (nontops.input.jet_pt.size() != nontops.input.jet_m.size()) {
        cout << "Unequal number of jet pts and jet pT widths!" << endl;
        cout << "there are " << nontops.input.jet_pt.size() << " jet pts and "
             << nontops.input.jet_m.size() << " jet masses" << endl;
        return;
    }
    if (nontops.input.jet_pt.size() != nontops.input.jet_pt_width.size()) {
        cout << "Unequal number of jet pts and jet pT widths!" << endl;
        cout << "there are " << nontops.input.jet_pt.size() << " jet pts and "
             << nontops.input.jet_pt_width.size() << " jet pt widths" << endl;
        return;
    }
    if (nontops.input.jet_pt.size() != nontops.input.jet_phi_width.size()) {
        cout << "Unequal number of jet pts and jet phi widths!" << endl;
        cout << "there are " << nontops.input.jet_pt.size() << " jet pts and "
             << nontops.input.jet_phi_width.size() << " jet phi widths" << endl;
        return;
    }
    if (nontops.input.jet_pt.size() != nontops.input.jet_eta_width.size()) {
        cout << "Unequal number of jet pts and jet eta widths!" << endl;
        cout << "there are " << nontops.input.jet_pt.size() << " jet pts and "
             << nontops.input.jet_eta_width.size() << " jet eta widths" << endl;
        return;
    }
    return;
}


inline void lightJetChiSquareMinimumSolver::Eval_covariance(
    recoil_minimizer_data &da)
{
    if (debug)
        cout << "lightJetChiSquareMinimumSolver::Eval_covariance\n";
    for (int i = 0; i < nontops.input.n_objs; i++) {
        if (debug) {
            cout << "Setting cartesian widths:" << endl;
            cout << "jetPtWidth = " << (nontops.input.jet_pt_width)[i] << endl;
            cout << "jetPhiWidth = " << (nontops.input.jet_phi_width)[i] << endl;
            cout << "jetEtaWidth = " << (nontops.input.jet_eta_width)[i] << endl;
        }
        const double pT = (nontops.input.jet_pt)[i];
        const double phi = (nontops.input.jet_phi)[i];
        const double eta = (nontops.input.jet_eta)[i];

        TMatrixD cov(3, 3);
        TMatrixD R(3, 3);
        TMatrixD RT(3, 3);

        cov[0][0] = pow((nontops.input.jet_pt_width)[i], 2);
        cov[0][1] = 0;
        cov[0][2] = 0;
        cov[1][0] = 0;
        cov[1][1] = pow((nontops.input.jet_phi_width)[i], 2);
        cov[1][2] = 0;
        cov[2][0] = 0;
        cov[2][1] = 0;
        cov[2][2] = pow((nontops.input.jet_eta_width)[i], 2);
        R[0][0] = cos(phi);
        R[0][1] = -pT * sin(phi);
        R[0][2] = 0;
        R[1][0] = sin(phi);
        R[1][1] = pT * cos(phi);
        R[1][2] = 0;
        //R[2][0] = tan(eta);
        R[2][0] = sinh(eta);
        R[2][1] = 0;
        //R[2][2] = pT * (1 + pow(tan(eta), 2));
        R[2][2] = pT * cosh(eta);
        RT.Transpose(R);

        // tranform to cartesian coordinate system
        da.cov[i] = R * cov * RT;
        da.cov_rad[i] = cov;

        if (debug) {
            cout << "Calculating widths:\n"
                 << "pT: " << pT << " +- " << (nontops.input.jet_pt_width)[i] << endl
                 << "phi: " << phi << " +- " << (nontops.input.jet_phi_width)[i] << endl
                 << "eta: " << eta << " +- " << (nontops.input.jet_eta_width)[i] << endl
                 << "px: " << pT * cos(phi) << " +- " << sqrt(cov[0][0]) << endl
                 << "py: " << pT * sin(phi) << " +- " << sqrt(cov[1][1]) << endl
                 << "pz: " << pT * tan(eta) << " +- " << sqrt(cov[1][1])
                 << endl;

            cout << "Covariance matrix in cartesian coords:\n";
            da.cov[i].Print();
        }
    }
}

inline void
lightJetChiSquareMinimumSolver::Eval_cov_sum(recoil_minimizer_data &da)
{
    da.inv_sum_cov.ResizeTo(3, 3);
    da.inv_sum_cov.Zero();
    for (int i = 0; i < nontops.input.n_objs; ++i)
        da.inv_sum_cov += da.cov[i];

    dynamic_cast<TDecompLU *>(inverter3D_)->SetMatrix(TMatrixD(da.inv_sum_cov));
    // checkDecomp = inverter3D_->Decompose();
    dynamic_cast<TDecompLU *>(inverter3D_)->Invert(da.inv_sum_cov);
    // da.inv_sum_cov3D.Print();

    for (int i = 0; i < nontops.input.n_objs; ++i)
        da.inv_sum_x_cov[i] = da.cov[i] * da.inv_sum_cov;
}

void lightJetChiSquareMinimumSolver::calcMin()
{
    if (do3D_) {
        if (data_.dxCheck == recoil_px(bigstruct) && data_.dyCheck == recoil_py(bigstruct) &&
            data_.dzCheck == recoil_pz(bigstruct))
            return;
    } else {
        if (data_.dxCheck == recoil_px(bigstruct) && data_.dyCheck == recoil_py(bigstruct))
            return;
    }

    // cout << "Calculating minimum chi^2" << endl;

    data_.dxCheck = recoil_px(bigstruct);
    data_.dyCheck = recoil_py(bigstruct);
    if (do3D_)
        data_.dzCheck = recoil_pz(bigstruct);
    else
        data_.dzCheck = 0;

    const double dp_arr[3] = {data_.dxCheck, data_.dyCheck, data_.dzCheck};
    const TVectorD dp(3, dp_arr);

    for (int i = 0; i < nontops.input.n_objs; ++i) {
        const TVectorD dp_i = data_.inv_sum_x_cov[i] * dp;

        nontops.best_innermost_params.jet_dif_px[i] = dp_i[0];
        nontops.best_innermost_params.jet_dif_py[i] = dp_i[1];
        nontops.best_innermost_params.jet_dif_pz[i] = dp_i[2];
    }

    double chi2 = dp * (data_.inv_sum_cov * dp);
    nontops.best_innermost_params.chi2 = chi2;
}

void lightJetChiSquareMinimumSolver::printResults()
{
    for (int i = 0; i < nontops.input.n_objs; ++i) {
        cout << "delta px " << i + 1 << " = " << nontops.best_innermost_params.jet_dif_px[i] << endl;
        cout << "delta py " << i + 1 << " = " << nontops.best_innermost_params.jet_dif_py[i] << endl;
        cout << "delta pz " << i + 1 << " = " << nontops.best_innermost_params.jet_dif_pz[i] << endl;
    }
}


#endif
