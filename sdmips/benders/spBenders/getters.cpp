#include "spbenders.h"

Master &spBenders::get_master(int stage, int node)
{
  d_master = make_unique<Master>(node_data(stage, outcome(stage, node)),
                                 stage == d_data.nstages() - 1,
                                 d_env);


  if (stage < d_cuts.size())
  {
    for (Cut const &cut : apx(stage, node))
      d_master->push_cut(cut);
  }

  return *d_master;
}

v_enum &spBenders::get_enums(int stage, int node)
{
  d_gens.clear();

  int future = stage + 1;
  d_gens.reserve(d_data.outcomes(future));

  vpath path(stage + 2);
  iota(path.begin(), path.end(), 0);

  vector<NodeData> nodes;
  nodes.reserve(path.size());
  int box = max(stage - d_depth + 1, 0);
  for (int lvl = 0; lvl < box; ++lvl)
  {
    nodes.push_back(node_data(lvl, 0));
    nodes[lvl].to_box();
  }

  vector<int> sub_path = tail(stage, node);
  for (int depth = 0; depth != sub_path.size(); ++depth)
    nodes.push_back(node_data(box + depth, sub_path[depth]));

  nodes.resize(nodes.size() + 1);

  bool leaf = (future == d_data.nstages() - 1);
  for (size_t out = 0; out != d_data.outcomes(future); ++out)
  {
    nodes.back() = node_data(future, out);
    d_gens.emplace_back(Enumerator{nodes, path, path.size() - 1, leaf, d_env});
  }

  for (Enumerator &gen : d_gens)
  {
    for (Cut const &cut : apx(stage, node))
      gen.add_cut(cut);
  }

  if (future == d_cuts.size())     // stage = T - 1: enumerators are leafs
    return d_gens;
    // add contemporary cuts
  vector<int> childs = children(stage, node);
  for (int idx = 0; idx != childs.size(); ++idx)
  {
    int child = childs[idx];
    for (Cut const &cut : apx(future, child))
      d_gens[idx].add_cut(cut);
  }

  return d_gens;
}

outer_apx const &spBenders::apx(int stage, int node)
{
  return d_depth == 0 ? d_cuts[stage].front() : d_cuts[stage][node];
}