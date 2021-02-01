#include "tree.h"

Tree::Tree(Stagewise const &sw)
{
  vector<int> parents { -1 };
  for (stage_data const &stage : sw.d_stages)
  {
    vector<int> children;
    size_t nparents = parents.size();
    for (int parent : parents)
    {
      for (NodeData data : stage)
      {
        data.d_prob /= nparents;
        children.push_back(add_node(data, parent));
      }
    }
    parents = children;
    children.clear();
  }
}
