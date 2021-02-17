#include "dbenders.h"

dBenders::dBenders(GRBEnv &env, Stagewise &data, int depth)
:
Benders(env, data, depth)
{
  bool stage_specific = (depth == 0);

  size_t nstages = d_data.nstages();
  d_nodes.reserve(nstages);

  vpath epath {0};                                     // for constructing the vertex enumerators
  vector<NodeData> edata {node_data(0,0)};

  for (int stage = 0; stage != nstages; ++stage)
  {
    epath.push_back(stage + 1);
    edata.resize(epath.size());

    int start = max(stage - depth + 1, 0);
    if (stage_specific) start = stage;

    vector<vpath> sub_paths = enumerate_paths(start, stage);

    vnode nodes;
    nodes.reserve(sub_paths.size());
    bool leaf = (stage == nstages - 1);
    bool preleaf = (stage == nstages - 2);

    for (vpath &tail : sub_paths)
    {
      Master master {node_data(stage, tail.back()), leaf, d_env};
      if (leaf)
      {
        nodes.emplace_back(node{move(master), nullptr});
        continue;
      }

      if (not stage_specific)
      {
        for (size_t lvl = 0; lvl != tail.size(); ++lvl)
          edata[start + lvl] = node_data(start + lvl, tail[lvl]);
      }

      int future = stage + 1;
      if (not stage_specific or nodes.empty())
      {
        shared_ptr<v_enum> e_ptr = make_shared<v_enum>();
        e_ptr->reserve(d_data.outcomes(future));

        for (int out = 0; out != d_data.outcomes(future); ++out)
        {
          edata.back() = node_data(future, out);
          e_ptr->emplace_back(Enumerator(edata, epath, epath.size() - 1, preleaf, d_env));
        }

        nodes.emplace_back(node{move(master), e_ptr});
      } else
        nodes.emplace_back(node{move(master), nodes.front().second});
    }
    d_nodes.emplace_back(move(nodes));

    for (NodeData &node : edata)
      node.to_box();
  }
}
