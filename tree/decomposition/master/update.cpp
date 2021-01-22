#include "master.h"

void Master::update(Solution const &sol)
{
  int ncons = d_data.ncons();
  int ncuts = d_cuts.size();

  arma::vec rhs_con = d_data.d_rhs - d_data.d_Bmat.t() * sol.d_x[0];

  vector<double> rhs(rhs_con.memptr(), rhs_con.memptr() + ncons);
  rhs.reserve(ncons + ncuts);

  for (Cut &cut : d_cuts)
    rhs.push_back(compute_rhs(cut, sol));

  GRBConstr *mip_cons = d_mip.getConstrs();
  d_mip.set(GRB_DoubleAttr_RHS, mip_cons, rhs.data(), rhs.size());
  delete[] mip_cons;

  GRBConstr *lp_cons = d_lp.getConstrs();
  d_lp.set(GRB_DoubleAttr_RHS, lp_cons, rhs.data(), rhs.size());
  delete[] lp_cons;

  d_mip.update();
  d_lp.update();
}