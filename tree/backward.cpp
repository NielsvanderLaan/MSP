#include "tree.h"

bool Tree::backward(int node)
{
  bool ret = false;

  if (is_leaf(node))
    return ret;

  for (int child : d_children[node])
  {
    if (backward(child))      // cut was added in one of the descendants
      ret = true;
  }

  //Cut cut = scaled_cut(node);     // TODO: choose by argument. Use abstract classes?
  //Cut cut = cpt_scaled_cut(node);
  Cut cut = sddp_cut(node);

  cout << "node: " << node << '\n';
  cut.print();
  d_masters[node].lp_forward().print();

  if (add_cut(node, cut))
    ret = true;




  return ret;
}