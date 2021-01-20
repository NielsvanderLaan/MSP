#include "lsde.h"

void LSDE::add(int id)
{
  NodeData data = d_tree.d_nodes[id];

  arma::vec costs = data.d_prob * data.d_costs;

  GRBVar *vars = d_model.addVars(data.d_lb.memptr(),
                                 data.d_ub.memptr(),
                                 costs.memptr(),
                                 data.d_types.memptr(),
                                 nullptr,
                                 data.d_costs.n_elem);

  d_vars[id] = vector<GRBVar>(vars, vars + data.d_costs.n_elem);

  GRBLinExpr lhs[data.ncons()];
  for (auto iter = data.d_Amat.begin(); iter != data.d_Amat.end(); ++iter)
    lhs[iter.col()] += *iter * vars[iter.row()];

  if (id > 0)
  {
    vector<GRBVar> prev_vars = d_vars[d_tree.d_ancestor[id]];

    for (auto iter = data.d_Bmat.begin(); iter != data.d_Bmat.end(); ++iter)
      lhs[iter.col()] += *iter * prev_vars[iter.row()];
  }

  GRBConstr *constrs = d_model.addConstrs(lhs,
                                          data.d_senses.memptr(),
                                          data.d_rhs.memptr(),
                                          nullptr,
                                          data.ncons());

  delete[] vars;
  delete[] constrs;

  d_model.update();
}