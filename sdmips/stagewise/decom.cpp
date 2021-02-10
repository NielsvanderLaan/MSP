#include "stagewise.h"

void Stagewise::decom(GRBEnv &env, int depth)
{
  d_depth = depth;
  if (depth == 0)
    return decom(env);

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
    vector<vpath> sub_paths = enumerate_paths(start, stage);

    vnode nodes;
    nodes.reserve(sub_paths.size());
    bool leaf = (stage == nstages - 1);
    bool preleaf = (stage == nstages - 2);
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

      v_enum *e_ptr = new v_enum;
      e_ptr->reserve(outcomes(stage + 1));
      for (int out = 0; out != outcomes(stage + 1); ++out)
      {
        edata.back() = d_stages[stage + 1][out];
        e_ptr->emplace_back(Enumerator(edata, epath, epath.size() - 1, preleaf, env));
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

void Stagewise::decom(GRBEnv &env)
{
  size_t nstages = d_stages.size();
  d_nodes.reserve(nstages);

  vpath fpath {};            // for constructing the Fenchel subproblems
  vector<NodeData> fdata;
  vpath epath {0};           // for constructing the vertex enumerators
  vector<NodeData> edata {d_stages[0][0]};

  for (int stage = 0; stage != nstages; ++stage)
  {
    fpath.push_back(stage);
    fdata.resize(fpath.size());
    epath.push_back(stage + 1);
    edata.resize(epath.size());

    vnode nodes;
    nodes.reserve(outcomes(stage));

    bool leaf = (stage == nstages - 1);
    v_enum *e_ptr = new v_enum;           // shared enumerator objects
    for (size_t out = 0; out != outcomes(stage); ++out)
    {
      NodeData const& data = d_stages[stage][out];
      Master master {data, leaf, env};
      if (leaf)
      {
        nodes.emplace_back(node{move(master), Enumerator{}, nullptr});
        continue;
      }

      if (out == 0)
      {
        int future = stage + 1;
        e_ptr->reserve(outcomes(future));
        for (int child = 0; child != outcomes(future); ++child)
        {
          edata.back() = d_stages[future][child];
          e_ptr->emplace_back(Enumerator(edata,
                                         epath,
                                         epath.size() - 1,
                                         stage + 1 == nstages - 1,
                                         env));
        }
      }

      fdata.back() = data;
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