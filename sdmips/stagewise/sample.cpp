#include "stagewise.h"

vector<vpath> Stagewise::sample(size_t nsamples)
{
  vector<vpath> ret(nsamples);
  for (vpath &path : ret)
    path.reserve(d_stages.size() - 1);   // final stage is not important

  for (auto it = d_stages.begin(); it != d_stages.end() - 1; ++it)
  {
    vector<double> prob = probs(*it);
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
  return enumerate_paths(0, d_stages.size() - 2);
}

