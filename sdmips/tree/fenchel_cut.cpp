#include "tree.h"

Cut Tree::fenchel_cut(int node, bool affine, double tol)
{
  return d_fenchel[node].feas_cut(d_masters[node].forward(), affine, tol);
}
