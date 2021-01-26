#include "../tree.h"

vvar Tree::add_to_lsde(int node, GRBModel &model, vvar const &parent_vars)
{
  NodeData &data = d_nodes[node];

  arma::vec costs = data.d_prob * data.d_costs;
  GRBVar *vars = model.addVars(data.d_lb.memptr(),
                                data.d_ub.memptr(),
                                costs.memptr(),
                                data.d_types.memptr(),
                                nullptr,
                                data.nvars());
  vvar xvars(vars, vars + data.nvars());
  delete[] vars;

  GRBLinExpr lhs[data.ncons()];
  for (auto iter = data.d_Amat.begin(); iter != data.d_Amat.end(); ++iter)
    lhs[iter.col()] += *iter * xvars[iter.row()];

  if (node > 0)
  {
    for (auto iter = data.d_Bmat.begin(); iter != data.d_Bmat.end(); ++iter)
      lhs[iter.col()] += *iter * parent_vars[iter.row()];
  }

  delete[] model.addConstrs(lhs,
                            data.d_senses.memptr(),
                            data.d_rhs.memptr(),
                            nullptr,
                            data.ncons());
  return xvars;
}
