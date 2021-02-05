#include "tree.h"

bool Tree::backward(bool affine, int node)
{
  bool ret = false;

  if (is_leaf(node))
    return ret;

  for (int child : d_children[node])
  {
    if (backward(affine,child))                    // cut was added in one of the descendants
      ret = true;
  }

  Cut cut = scaled_cut(node, affine);       // TODO: choose by argument. Use abstract classes?
  //Cut cut = cpt_scaled_cut(false, node);
  //Cut cut = sddp_cut(node);

  if (add_cut(node, cut))
    ret = true;

  return ret;
}