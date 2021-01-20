#include "master.h"

void Master::update(Solution &sol)
{
  int ncons = d_data.ncons();
  int ncuts = d_cuts.size();

  arma::vec rhs_con = d_data.d_rhs - d_data.d_Bmat.t() * sol.d_x[0];

  vector<double> rhs(rhs_con.memptr(), rhs_con.memptr() + ncons);
  rhs.reserve(ncons + ncuts);

  for (Cut &cut : d_cuts)
    rhs.push_back(compute_rhs(cut, sol));

  GRBConstr *cons = d_model.getConstrs();
  d_model.set(GRB_DoubleAttr_RHS, cons, rhs.data(), rhs.size());
  delete[] cons;

  d_model.update();
}