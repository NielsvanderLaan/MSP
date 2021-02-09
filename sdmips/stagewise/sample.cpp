#include "stagewise.h"

vector<vpath> Stagewise::sample(size_t nsamples)
{
  vector<vpath> ret(nsamples);

  for (stage_data const &stage : d_stages)
  {
    vector<double> prob = probs(stage);
    discrete_distribution<int> uni(prob.begin(), prob.end());
    for (vpath &path : ret)
      path.push_back(uni(d_engine));
  }

  return ret;
}

vector<vpath> Stagewise::enumerate_paths(int start, int end, vector<vpath> paths)
{
  vector<vpath> ret;
  int outcomes = d_stages[start].size();
  ret.reserve(paths.size() * outcomes);

  for (vpath &path : paths)
  {
    for (int outcome = 0; outcome != outcomes; ++outcome)
    {
      vpath copy = path;
      copy.push_back(outcome);
      ret.push_back(copy);
    }
  }

  if (start < end)
    return enumerate_paths(start + 1, end, ret);

  return ret;
}

vector<vpath> Stagewise::enumerate_paths()
{
  return enumerate_paths(0, d_stages.size() - 1);
}

