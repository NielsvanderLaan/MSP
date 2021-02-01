#include "tree.h"

Cut Tree::scaled_cut(int node, double tol)
{
  Cut ret;
  vector<int> path_nvars = nvars(node);

  double rho = d_masters[node].theta_n();
  double crho;

  init_enums(node);
  double par_prob = d_nodes[node].d_prob;
  do
  {
    ret = Cut(path_nvars);
    crho = -rho;

    for (int child : d_children[node])
    {
      double qnm = d_nodes[child].d_prob / par_prob;
      ret += qnm * d_enumerators[child].opt_cut(rho, tol);
      crho -= qnm * d_enumerators[child].crho();
    }
    rho += crho / (1 + ret.d_tau.back());
  } while (crho > tol);

  return ret;
}