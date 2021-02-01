#include "nodedata.h"

vvar NodeData::add_to_lsde(GRBModel &lsde, vvar const &parent_vars, double corr) const
{
  arma::vec costs = d_prob * d_costs * corr;
  GRBVar *vars = lsde.addVars(d_lb.memptr(),
                              d_ub.memptr(),
                              costs.memptr(),
                              d_types.memptr(),
                              nullptr,
                              nvars());
  vvar xvars(vars, vars + nvars());
  delete[] vars;

  GRBLinExpr lhs[ncons()];
  for (auto iter = d_Amat.begin(); iter != d_Amat.end(); ++iter)
    lhs[iter.col()] += *iter * xvars[iter.row()];

  if (parent_vars.size() > 0)
  {
    for (auto iter = d_Bmat.begin(); iter != d_Bmat.end(); ++iter)
      lhs[iter.col()] += *iter * parent_vars[iter.row()];
  }

  delete[] lsde.addConstrs(lhs,
                           d_senses.memptr(),
                           d_rhs.memptr(),
                           nullptr,
                           ncons());
  return xvars;

}