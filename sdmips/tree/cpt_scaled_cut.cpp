#include "tree.h"

Cut Tree::cpt_scaled_cut(int node, bool affine, double tol)
{
  Cut ret;
  vector<int> path_nvars = nvars(node);

  double rho = d_masters[node].theta_n();
  double crho;

  double par_prob = d_nodes[node].d_prob;

  do
  {
    ret = Cut(path_nvars);
    crho = -rho;

    for (int child : d_children[node])
    {
      Master &master = d_masters[child];
      double qnm = d_nodes[child].d_prob / par_prob;

      master.set_rho(rho);
      solve_master(child, affine, false, false);      // SND MIP using Fenchel cuts
      ret += qnm * master.opt_cut();                           // LP duality --> optimality cut
      crho -= qnm * master.lp_obj();
    }

    rho += crho / (1 + ret.d_tau.back());
  } while (crho > tol);

  return ret;
}