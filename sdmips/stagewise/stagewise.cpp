#include "stagewise.h"

void Stagewise::add_node(NodeData const &data)
{
  int stage = data.d_stage;
  if (d_stages.size() < stage)
    d_stages.resize(stage);

  if (stage == 1)
    assert(d_stages[0].empty());

  d_stages[stage - 1].push_back(data);
}

GRBModel Stagewise::lsde(GRBEnv &env)
{
  GRBModel model(env);
  vector<vvar> parent_vars { vvar{} };
  vector<vvar> children_vars;

  for (stage_data const &stage : d_stages)
  {
    double corr = 1.0 / parent_vars.size();
    for (NodeData const& data : stage)
    {
      for (vvar const &parent : parent_vars)
        children_vars.emplace_back(data.add_to_lsde(model, parent, corr));
    }
    parent_vars = children_vars;
    children_vars.clear();
  }

  return model;
}

Stagewise::~Stagewise()
{
  if (d_nodes.empty())
    return;

  for (auto it = d_nodes.begin(); it != d_nodes.end() - 1; ++it)
  {
    if (d_depth > 0)
    {
      for (node &sub: *it)
        delete get<2>(sub);
    }
    else                                  // enumerators are shared
      delete get<2>(it->front());
  }
}

