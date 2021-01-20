#include "tree.h"

Cut Tree::lp_cut(int id, Solution &sol)
{
  Cut ret(nvars(id));
  double par_prob = d_nodes[id].d_prob;

  for (int child : d_children[id])
    ret += d_masters[child].compute_cut(sol) * (d_nodes[child].d_prob / par_prob);

  return ret;
}