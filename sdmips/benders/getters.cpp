#include "benders.h"

NodeData &Benders::node_data(int stage, int outcome)
{
  return d_data.d_stages[stage][outcome];
}

int Benders::master_idx(int stage, vpath const &path) const
{
  int ret = path[stage];
  int jump = 1;
  for (int lvl = 1; lvl < min(d_depth, stage + 1); ++lvl)
  {
    jump *= d_data.outcomes(stage - lvl + 1);
    ret += path[stage - lvl] * jump;
  }

  return ret;
}

vector<int> Benders::parents(int stage, int node) const
{
  assert(d_depth > 0);

  int npaths = d_data.outcomes(max(stage - d_depth, 0));
  vector<int> ret(npaths);

  int offset = 0;
  int jump = 1;
  vpath sub_path = tail(stage, node);

  for (int lvl = 1; lvl < sub_path.size(); ++lvl)
  {
    sub_path.pop_back();
    offset += sub_path.back() * jump;
    jump *= d_data.outcomes(stage- lvl);
  }

  for (int idx = 0; idx != npaths; ++idx)
    ret[idx] = idx * jump + offset;

  return ret;
}

vector<int> Benders::children(int stage, int node) const
{
  int n_outcomes = d_data.outcomes(stage + 1);
  vector<int> ret(n_outcomes);
  if (stage < d_depth)
  {
    iota(ret.begin(), ret.end(), node * n_outcomes);
    return ret;
  }

  int skip = 1;
  for (int depth = 0; depth < d_depth - 1; ++depth)
    skip *= d_data.outcomes(stage - depth);

  iota(ret.begin(), ret.end(), (node % skip) * n_outcomes);
  return ret;
}

vector<int> Benders::tail(int stage, int node) const
{
  vector<int> ret(min(stage + 1, d_depth));

  for (size_t depth = 0; depth != ret.size(); ++depth)
  {
    ret[depth] = outcome(stage - depth, node);
    node /= d_data.outcomes(stage - depth);
  }
  reverse(ret.begin(), ret.end());
  return ret;
}

int Benders::outcome(int stage, int node) const
{
  return node % d_data.outcomes(stage);
}