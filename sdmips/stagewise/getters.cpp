#include "stagewise.h"
#include "stagewise.h"

int Stagewise::nstages() const
{
  return d_stages.size();
}

int Stagewise::outcomes(int stage) const
{
  return d_stages[stage].size();
}

vector<int> Stagewise::nvars(int stage) const
{
  vector<int> ret(stage + 1);
  for (int lvl = 0; lvl != stage + 1; ++lvl)
    ret[lvl] = d_stages[lvl].front().nvars();
  return ret;
}

vector<double> Stagewise::probs(int stage) const
{
  vector<double> probs;
  probs.reserve(outcomes(stage));
  for (NodeData const &data : d_stages[stage])
    probs.push_back(data.d_prob);
  return probs;
}

