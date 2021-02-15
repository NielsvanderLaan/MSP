#include "master.h"

bool Master::add_cut(Cut const &cut, double tol)
{
  double lhs_val = 0;
  double rhs_val = 0;

  Solution sol = forward();
  if (cut.depth() == sol.depth())
  {
    lhs_val = compute_lhs(cut, sol);
    rhs_val = compute_rhs(cut, sol);
  }

  if (lhs_val >= rhs_val - tol)
    return false;

  push_cut(cut, rhs_val);

  return true;
}

void Master::add(Cut const &cut)
{
  add_cut(cut, -GRB_INFINITY);
}

void Master::push_cut(Cut const &cut, double rhs)
{
  d_cuts.push_back(cut);

  GRBLinExpr lhs = (1 + cut.d_tau.back()) * d_lp_theta;
  lhs.addTerms(cut.d_beta.back().memptr(), d_lp_xvars.data(), d_lp_xvars.size());
  d_lp->addConstr(lhs >= rhs);
  d_lp->update();

  if (not cut.d_feas)
  {
    GRBLinExpr mip_lhs = (1 + cut.d_tau.back()) * d_mip_theta;
    mip_lhs.addTerms(cut.d_beta.back().memptr(), d_mip_xvars.data(), d_mip_xvars.size());
    d_mip->addConstr(mip_lhs >= rhs);
    d_mip->update();
  }
}

void Master::push_cut(Cut const &cut)
{
  push_cut(cut, cut.d_alpha);
}

