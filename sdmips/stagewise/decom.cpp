#include "stagewise.h"

void Stagewise::decom(GRBEnv &env, int depth)
{
  d_depth = depth;
  size_t nstages = d_stages.size();
  d_nodes.reserve(nstages);

  vpath fpath {};            // for constructing the Fenchel subproblems
  vector<NodeData> fdata;
  vpath epath {0};           // for constructing the vertex enumerators
  vector<NodeData> edata;

  for (int stage = 0; stage != nstages; ++stage)
  {
    fpath.push_back(stage);
    fdata.resize(fpath.size());
    epath.push_back(stage + 1);
    edata.resize(epath.size());

    int start = max(stage - depth + 1, 0);
    if (depth == 0)
      start = stage;
    vector<vpath> sub_paths = enumerate_paths(start, depth);

    vnode nodes;
    nodes.reserve(sub_paths.size());
    bool leaf = (stage == nstages - 1);
    for (vpath &tail : sub_paths)
    {
      Master master {d_stages[stage][tail.back()], leaf, env};
      if (leaf)
      {
        nodes.emplace_back(node{move(master), Enumerator{}, nullptr});
        continue;
      }

      for (size_t lvl = 0; lvl != tail.size(); ++lvl)
      {
        int st = start + lvl;
        int out = tail[lvl];
        fdata[st] = d_stages[st][out];
        edata[st] = d_stages[st][out];
      }

      v_enum *e_ptr;
      if (depth == 0 && nodes.size() > 0)
        e_ptr = &get_enums(stage, 0);
      else
      {
        e_ptr = new v_enum;
        for (int out = 0; out != outcomes(stage + 1); ++out)
        {
          edata.back() = d_stages[stage + 1][out];
          e_ptr->emplace_back(Enumerator(edata, epath, epath.size() - 1, leaf, env));
        }
      }

      nodes.emplace_back(node{move(master),
                              Enumerator {fdata, fpath, fpath.size(), leaf, env},
                              e_ptr});
    }
    d_nodes.emplace_back(move(nodes));

    for (NodeData &data : fdata)
      data.to_box();
    for (NodeData &data : edata)
      data.to_box();
  }
}