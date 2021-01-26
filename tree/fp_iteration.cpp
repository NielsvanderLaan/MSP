#include "tree.h"

Cut Tree::fp_iteration(int node, vector<CutFamily> &gens, double tol)
{
  Cut ret;
  vector<int> path_nvars = nvars(node);

  if (is_leaf(node))
    return Cut{};

  /*
  for (CutFamily &gen : gens)
    gen.set();
  */

  int child = d_children[node][0];
  double rho = d_masters[child].d_state.d_theta.back();

  double crho;

  double par_prob = d_nodes[node].d_prob;

  do
  {
    ret = Cut(path_nvars);
    crho = -rho;

    for (int child : d_children[node])
    {
      double qnm = d_nodes[child].d_prob / par_prob;
      ret += qnm * gens[child].compute_cut(rho, tol);
      crho -= qnm * gens[child].crho();
    }

    rho += crho / (1 + ret.d_tau.back());
  } while (crho > tol);

  return ret;
}







