#include "../tree/tree.h"

bool Tree::is_leaf(int node)
{
  return d_children[node].size() == 0;
}