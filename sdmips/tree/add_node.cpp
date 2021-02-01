#include "tree.h"

int Tree::add_node(NodeData const &node, int ancestor)
{
  int id = d_nodes.size();
  assert(ancestor < id);

  d_nodes.push_back(node);

  if (ancestor != -1)
    d_children[ancestor].push_back(id);
  else
    assert(id == 0);

  d_ancestor[id] = ancestor;

  return id;
}