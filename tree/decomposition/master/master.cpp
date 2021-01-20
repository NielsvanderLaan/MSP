#include "master.h"

Master::Master(NodeData &data, bool leaf, GRBEnv &env)
:
d_data(data),
d_model(d_data.to_model(env, true))
{
  GRBVar *vars = d_model.getVars();
  d_xvars = vector<GRBVar>(vars, vars + d_data.nvars());
  delete[] vars;

  if (leaf)
    return;

  d_theta = d_model.addVar(d_data.d_L, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
  d_model.update();
}

Master::Master(const Master &other)
:
d_data(other.d_data),
d_model(other.d_model),
d_cuts(other.d_cuts)
{
  GRBVar *vars = d_model.getVars();
  d_xvars = vector<GRBVar>(vars, vars + d_data.nvars());
  d_theta = vars[d_data.nvars()];

  delete[] vars;
}