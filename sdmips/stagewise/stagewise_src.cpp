#include "stagewise.h"

void Stagewise::decom(GRBEnv &env)
{
  d_masters.reserve(d_stages.size());

  for (stage_data const &stage : d_stages)
  {
    vmaster masters;
    for (NodeData const &data : stage)
      masters.emplace_back(Master{data, &stage == &d_stages.back(), env});

    d_masters.push_back(masters);
  }
}

void Stagewise::sddmip()
{
  size_t max_iter = 10;
  for (int iter = 0; iter != max_iter; ++iter)    // TODO: stopping criterion
  {
    vector<path> paths = sample();
    vector<vsol> sols = forward(paths);
    cout << d_masters[0][0].obj() << '\n';

    backward(sols);
  }
}

vector<path> Stagewise::sample(size_t nsamples)   // TODO
{
  size_t npaths = 1;
  for (auto &stage : d_stages)
    npaths *= stage.size();
  vector<path> paths(npaths);

  for (auto &stage : d_stages)
  {
    for (size_t path = 0; path != npaths; ++path)
      paths[path].push_back(path % stage.size());
  }

  return paths;
}

vector<vsol> Stagewise::forward(vector<path> &paths)
{
  size_t nstages = d_stages.size();
  vector<vsol> ret(nstages - 1);

  for (auto &path : paths)
  {
    for (int stage = 0; stage != nstages - 1; ++stage)
    {
      int node = path[stage];
      int child = path[stage + 1];

      solve(stage, node);
      Solution forward = d_masters[stage][node].forward();
      ret[stage].push_back(forward);

      d_masters[stage + 1][child].update(forward);
    }
  }
  ret[0].resize(1);
  return ret;
}

void Stagewise::solve(int stage, int node)
{
  d_masters[stage][node].solve_lp();
}

void Stagewise::backward(vector<vsol> const &sols)
{
  int stage = d_stages.size() - 1;

  for (auto it = sols.rbegin(); it != sols.rend(); ++it)
  {
    --stage;
    for (Solution const &sol : *it)
    {
      Cut cut = sddp_cut(stage, sol);     // add cuts in bulk? (scaled cuts depend on previous cuts)
      add_cut(cut, stage);
    }
  }
}

Cut Stagewise::sddp_cut(int stage, Solution const &sol)
{
  Cut ret(nvars(stage));
  vmaster subs = d_masters[stage + 1];

  for (int child = 0; child != subs.size(); ++child)
  {
    subs[child].update(sol);
    Cut cut = subs[child].opt_cut();
    cut.print();
    cout << d_stages[stage][child].d_prob << '\n';

    ret += d_stages[stage][child].d_prob *  cut;
  }
  return ret;
}

void Stagewise::add_cut(Cut &cut, int stage)
{
  for (Master &master : d_masters[stage])
    master.add(cut);

  // TODO: update enumerators etc.
}



GRBModel Stagewise::lsde(GRBEnv &env)
{
  GRBModel model(env);
  vector<vvar> parent_vars { vvar{} };
  vector<vvar> children_vars;

  for (stage_data const &stage : d_stages)
  {
    double corr = 1.0 / parent_vars.size();
    for (NodeData const& data : stage)
    {
      for (vvar const &parent : parent_vars)
        children_vars.emplace_back(data.add_to_lsde(model, parent, corr));
    }
    parent_vars = children_vars;
    children_vars.clear();
  }

  return model;
}

void Stagewise::add_node(NodeData const &data)
{
  int stage = data.d_stage;
  if (d_stages.size() < stage)
    d_stages.resize(stage);

  if (stage == 1)
    assert(d_stages[0].empty());

  d_stages[data.d_stage - 1].push_back(data);
}

vector<int> Stagewise::nvars(int stage) const
{
  vector<int> ret(stage + 1);
  for (int lvl = 0; lvl != stage + 1; ++lvl)
    ret[lvl] = d_stages[lvl][0].nvars();
  return ret;
}