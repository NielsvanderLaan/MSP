#include "master.h"

void Master::update(Solution const &sol)
{
  d_state = sol;

  int ncons = d_data.ncons();
  int ncuts = d_cuts.size();

  arma::vec rhs_con = d_data.d_rhs - d_data.d_Bmat.t() * sol.d_x.back();

  vector<double> lp_rhs(rhs_con.memptr(), rhs_con.memptr() + ncons);
  vector<double> mip_rhs = lp_rhs;
  lp_rhs.reserve(ncons + ncuts);

  for (Cut &cut : d_cuts)
  {
    double rhs_val = compute_rhs(cut, sol);
    lp_rhs.push_back(rhs_val);

    if (not cut.d_feas)
      mip_rhs.push_back(rhs_val);
  }

  GRBConstr *mip_cons = d_mip.getConstrs();
  d_mip.set(GRB_DoubleAttr_RHS, mip_cons, mip_rhs.data(), mip_rhs.size());
  delete[] mip_cons;

  GRBConstr *lp_cons = d_lp.getConstrs();
  d_lp.set(GRB_DoubleAttr_RHS, lp_cons, lp_rhs.data(), lp_rhs.size());
  delete[] lp_cons;

  d_mip.update();
  d_lp.update();
}