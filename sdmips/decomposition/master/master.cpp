#include "master.h"

Master::Master(NodeData const &data, bool leaf, GRBEnv &env)
:
d_data(data),
d_mip(d_data.to_model(env, false)),
d_lp(d_data.to_model(env, true))
{
  GRBVar *mip_vars = d_mip.getVars();
  d_mip_xvars = vvar(mip_vars, mip_vars + d_data.nvars());
  delete[] mip_vars;

  GRBVar *lp_vars = d_lp.getVars();
  d_lp_xvars = vvar(lp_vars, lp_vars + d_data.nvars());
  delete[] lp_vars;

      /*
       * lb and ub on theta (representing expected costs in future stages)
       * leaf node --> theta = 0
       */
  double lb = leaf ? 0.0 : d_data.d_L;
  double ub = leaf ? 0.0 : GRB_INFINITY;

  d_mip_theta = d_mip.addVar(lb, ub, 1.0, GRB_CONTINUOUS);
  d_lp_theta = d_lp.addVar(lb, ub, 1.0, GRB_CONTINUOUS);

  d_mip.update();
  d_lp.update();
}

Master::Master(const Master &other)
:
  d_data(other.d_data),
  d_mip(other.d_mip),
  d_lp(other.d_lp),
  d_cuts(other.d_cuts),
  d_state(other.d_state),
  d_x_n(other.d_x_n),
  d_theta_n(other.d_theta_n)
{
  int nvars = d_data.nvars();

  GRBVar *mip_vars = d_mip.getVars();
  d_mip_xvars = vector<GRBVar>(mip_vars, mip_vars + nvars);
  d_mip_theta = mip_vars[nvars];
  delete[] mip_vars;

  GRBVar *lp_vars = d_lp.getVars();
  d_lp_xvars = vvar(lp_vars, lp_vars + nvars);
  d_lp_theta = lp_vars[nvars];
  delete[] lp_vars;
}