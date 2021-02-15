#include "stagewise.h"

void Stagewise::decom(int depth)
{
  bool stage_specific = (depth == 0);
  d_depth = depth;

  size_t nstages = d_stages.size();
  d_nodes.reserve(nstages);

  vpath epath {0};                            // for constructing the vertex enumerators
  vector<NodeData> edata {d_stages[0][0]};

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
      Master master {d_stages[stage][tail.back()], leaf, d_env};
      if (leaf || true)
      {
        nodes.emplace_back(node{move(master), nullptr});
        continue;
      }

      if (not stage_specific)
      {
        for (size_t lvl = 0; lvl != tail.size(); ++lvl)
          edata[start + lvl] = d_stages[start + lvl][tail[lvl]];
      }

      if (not stage_specific or nodes.empty())
      {
        v_enum *e_ptr = new v_enum;
        e_ptr->reserve(outcomes(stage + 1));
        for (int out = 0; out != outcomes(stage + 1); ++out)
        {
          edata.back() = d_stages[stage + 1][out];
          e_ptr->emplace_back(Enumerator(edata, epath, epath.size() - 1, preleaf, d_env));
        }
        nodes.emplace_back(node{move(master), e_ptr});
      } else
      {
        nodes.emplace_back(node{move(master), nodes.front().second});
      }

    }
    d_nodes.emplace_back(move(nodes));

    for (NodeData &data : edata)
      data.to_box();
  }
}