#include "tree.h"

Cut Tree::lp_cut(int id, Solution const &sol)
{
  Cut ret(nvars(id));
  double par_prob = d_nodes[id].d_prob;

  for (int child : d_children[id])
    ret += (d_nodes[child].d_prob / par_prob) * d_masters[child].compute_cut(sol);

  return ret;
}