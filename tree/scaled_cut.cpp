#include "tree.h"

Cut Tree::scaled_cut(int node, Solution sol, double tol)
{
  Cut ret;
  vector<int> path_nvars = nvars(node);

  double rho = sol.d_theta.back();
  double crho = tol + 1;

  init_enums(node, sol);
  double par_prob = d_nodes[node].d_prob;

  while (crho > tol)
  {
    ret = Cut(path_nvars);

    crho = -rho;

    for (int child : d_children[node])
    {
      double qnm = d_nodes[child].d_prob / par_prob;
      ret += qnm * d_enumerators[child].generate_cut(rho);
      crho -= qnm * d_enumerators[child].crho();
    }

    rho += crho / (1 + ret.d_tau.back());
  }

  return ret;
}