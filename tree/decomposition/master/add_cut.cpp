#include "master.h"

bool Master::add_cut(Cut &cut, Solution &sol, double tol)
{
  double lhs_val = compute_lhs(cut, sol);
  double rhs_val = compute_rhs(cut, sol);

  if (lhs_val >= rhs_val - tol)
    return false;

  d_cuts.push_back(cut);

  GRBLinExpr lhs = (1 + cut.d_tau[0]) * d_theta;
  lhs.addTerms(cut.d_beta[0].memptr(), d_xvars.data(), d_xvars.size());

  d_model.addConstr(lhs >= rhs_val);

  d_model.update();

  return true;
}