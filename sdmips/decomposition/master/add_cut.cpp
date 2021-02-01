#include "master.h"

bool Master::add_cut(Cut const &cut, double tol)
{
  Solution sol = forward();

  double lhs_val = compute_lhs(cut, sol);
  double rhs_val = compute_rhs(cut, sol);

  if (lhs_val >= rhs_val - tol)
    return false;

  d_cuts.push_back(cut);

  GRBLinExpr lhs = (1 + cut.d_tau.back()) * d_lp_theta;
  lhs.addTerms(cut.d_beta.back().memptr(), d_lp_xvars.data(), d_lp_xvars.size());
  d_lp.addConstr(lhs >= rhs_val);
  d_lp.update();

  if (not cut.d_feas)
  {
    GRBLinExpr mip_lhs = (1 + cut.d_tau.back()) * d_mip_theta;
    mip_lhs.addTerms(cut.d_beta.back().memptr(), d_mip_xvars.data(), d_mip_xvars.size());
    d_mip.addConstr(mip_lhs >= rhs_val);
    d_mip.update();
  }

  return true;
}

void Master::add(const Cut &cut)
{
  add_cut(cut, -GRB_INFINITY);
}