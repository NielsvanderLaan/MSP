#include "stagewise.h"

int Stagewise::master_idx(int stage, vpath const &path)
{
  int ret = path[stage];
  int jump = 1;
  for (int lvl = 1; lvl < min(d_depth, stage + 1); ++lvl)
  {
    jump *= outcomes(stage - lvl + 1);
    ret += path[stage - lvl] * jump;
  }

  return ret;
}


vector<int> Stagewise::parents(int stage, vector<int> const &path)
{
  assert(d_depth > 0);

  int offset = 0;
  int jump = 1;

  for (int lvl = 1; lvl < min(d_depth, stage + 1); ++lvl)
  {
    offset += path[stage - lvl] * jump;
    jump *= outcomes(stage- lvl);
  }

  int npaths = outcomes(max(stage - d_depth, 0));
  vector<int> ret(npaths);
  for (int path = 0; path != npaths; ++path)
    ret[path] = path * jump + offset;

  return ret;
}

vector<int> Stagewise::children(int stage, int node)
{
  int n_outcomes = outcomes(stage + 1);
  vector<int> ret(n_outcomes);
  if (stage < d_depth)
  {
    iota(ret.begin(), ret.end(), node * n_outcomes);
    return ret;
  }

  int skip = 1;
  for (int depth = 0; depth != d_depth - 1; ++depth)
    skip *= outcomes(stage - depth);

  iota(ret.begin(), ret.end(), (node % skip) * n_outcomes);
  return ret;
}

vector<int> Stagewise::nvars(int stage) const
{
  vector<int> ret(stage + 1);
  for (int lvl = 0; lvl != stage + 1; ++lvl)
    ret[lvl] = d_stages[lvl].front().nvars();
  return ret;
}

vector<double> Stagewise::probs(stage_data const &stage) const
{
  vector<double> probs;
  probs.reserve(stage.size() - 1);
  for (NodeData const &data : stage)
    probs.push_back(data.d_prob);
  return probs;
}

Master &Stagewise::get_master(int stage, int node)
{
  return get<0>(d_nodes[stage][node]);
}

Enumerator &Stagewise::get_fenchel(int stage, int node)
{
  return get<1>(d_nodes[stage][node]);
}

v_enum &Stagewise::get_enums(int stage, int node)
{
  return *get<2>(d_nodes[stage][node]);
}

Solution Stagewise::solution(int stage, int node)
{
  return get_master(stage, node).forward();
}

int Stagewise::outcomes(int stage) const
{
  return d_stages[stage].size();
}