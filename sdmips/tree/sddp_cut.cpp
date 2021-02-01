#include "tree.h"

Cut Tree::sddp_cut(int node)
{
  Cut ret(nvars(node));
  double par_prob = d_nodes[node].d_prob;

  for (int child : d_children[node])
    ret += (d_nodes[child].d_prob / par_prob) * d_masters[child].opt_cut();

  return ret;
}