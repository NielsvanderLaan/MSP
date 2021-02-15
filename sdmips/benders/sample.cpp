#include "benders.h"

vector<vpath> Benders::sample(size_t nsamples)
{
  vector<vpath> ret(nsamples);
  for (vpath &path : ret)
    path.reserve(d_data.nstages() - 1);   // final stage is not important

  for (size_t stage = 0; stage != d_data.nstages() - 1; ++stage)
  {
    vector<double> prob = d_data.probs(stage);
    discrete_distribution<int> uni(prob.begin(), prob.end());
    for (vpath &path : ret)
      path.push_back(uni(d_engine));
  }

  return ret;
}

vector<vpath> Benders::enumerate_paths(int start, int end, vector<vpath> const &paths)
{
  vector<vpath> ret;
  int outcomes = d_data.outcomes(start);
  ret.reserve(paths.size() * outcomes);

  for (vpath const &path : paths)
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

vector<vpath> Benders::enumerate_paths()
{
  return enumerate_paths(0, d_data.nstages() - 2);
}

