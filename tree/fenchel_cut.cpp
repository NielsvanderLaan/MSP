#include "tree.h"

Cut Tree::fenchel_cut(int node, double tol)
{
  return d_fenchel[node].feas_cut(d_masters[node].lp_forward(), tol);
}
