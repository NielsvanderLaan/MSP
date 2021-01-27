#include "nodedata.h"

GRBModel NodeData::to_model(GRBEnv &env, bool lp)
{
  GRBModel model(env);

  GRBVar *vars = model.addVars(d_lb.memptr(),
                               d_ub.memptr(),
                               d_costs.memptr(),
                               lp ? nullptr : d_types.memptr(),
                               nullptr,
                                nvars());

  GRBLinExpr lhs[d_Amat.n_cols];
  for (auto iter = d_Amat.begin(); iter != d_Amat.end(); ++iter)
    lhs[iter.col()] += *iter * vars[iter.row()];

  GRBConstr *constrs = model.addConstrs(lhs,
                                        d_senses.memptr(),
                                        d_rhs.memptr(),
                                        nullptr,
                                        ncons());

  delete[] vars;
  delete[] constrs;

  model.update();

  return model;
}