#include "master.h"

bool Master::add_cut(Cut &cut, Solution &sol, double tol)
{
  double lhs_val = compute_lhs(cut, sol);
  double rhs_val = compute_rhs(cut, sol);

  if (lhs_val >= rhs_val - tol)
    return false;

  d_cuts.push_back(cut);

  GRBLinExpr mip_lhs = (1 + cut.d_tau[0]) * d_mip_theta;
  mip_lhs.addTerms(cut.d_beta[0].memptr(), d_mip_xvars.data(), d_mip_xvars.size());
  d_mip.addConstr(mip_lhs >= rhs_val);

  GRBLinExpr lhs = (1 + cut.d_tau[0]) * d_lp_theta;
  lhs.addTerms(cut.d_beta[0].memptr(), d_lp_xvars.data(), d_lp_xvars.size());
  d_lp.addConstr(lhs >= rhs_val);

  d_mip.update();
  d_lp.update();

  return true;
}