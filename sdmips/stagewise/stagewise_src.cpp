#include "stagewise.h"

void Stagewise::decom(GRBEnv &env)
{
  size_t nstages = d_stages.size();
  d_masters.reserve(nstages);
  d_enumerators.reserve(nstages);
  d_fenchel.reserve(nstages);

  vpath path {};
  vector<NodeData> nodes;
  for (stage_data const &stage : d_stages)
  {
    path.push_back(path.size());
    nodes.resize(path.size());

    vmaster masters;
    v_enum enums;
    v_enum fenchels;
    size_t outcomes = stage.size();
    masters.reserve(outcomes);
    enums.reserve(outcomes);
    fenchels.reserve(outcomes);

    bool leaf = (&stage == &d_stages.back());
    for (NodeData const &data : stage)
    {
      nodes.back() = data;

      masters.emplace_back(Master {data, leaf, env});
      enums.emplace_back(Enumerator {nodes, path, path.size() - 1, leaf, env});
      fenchels.emplace_back(Enumerator {nodes, path, path.size(), leaf, env});
    }

    d_masters.emplace_back(move(masters));
    d_enumerators.emplace_back(move(enums));
    d_fenchel.emplace_back(move(fenchels));

    nodes.back().to_box();
  }
}

void Stagewise::sddmip(bool affine)
{
  size_t max_iter = 25;
  for (int iter = 0; iter != max_iter; ++iter)          // TODO: stopping criterion
  {
    vector<vpath> paths = sample(5);
    vector<vsol> sols = forward(paths, affine, false);
    cout << d_masters[0][0].obj() << endl;

    backward(sols, affine);
  }
}

vector<vsol> Stagewise::forward(vector<vpath> &paths, bool affine, bool lp)
{
  size_t nstages = d_stages.size();
  vector<vsol> ret(nstages - 1);

  for (auto &path : paths)
  {
    for (int stage = 0; stage != nstages - 1; ++stage)
    {
      int node = path[stage];
      int child = path[stage + 1];

      solve(stage, node, affine, lp, true);
      Solution forward = d_masters[stage][node].forward();
      ret[stage].push_back(forward);

      d_masters[stage + 1][child].update(forward);
    }
  }
  ret.begin()->resize(1);
  return ret;
}

void Stagewise::backward(vector<vsol> const &sols, bool affine)
{
  int stage = d_stages.size() - 1;

  for (auto it = sols.rbegin(); it != sols.rend(); ++it)
  {
    --stage;
    vector<Cut> cuts;
    cuts.reserve(it->size());
    for (Solution const &sol : *it)
    {
      //cuts.push_back(sddp_cut(stage, sol));
      cuts.push_back(scaled_cut(stage, sol, affine));
    }

    for (Cut &cut : cuts)
      add_cut(cut, stage);
  }
}

Cut Stagewise::scaled_cut(int stage, const Solution &sol, bool affine, double tol)
{
  Cut ret;
  vector<int> path_nvars = nvars(stage);

  double rho = sol.theta_n();
  double crho;

  init_enums(stage + 1, sol);

  do
  {
    ret = Cut(path_nvars);
    crho = -rho;
    for (Enumerator &gen : d_enumerators[stage + 1])
    {
      ret += gen.d_data.d_prob * gen.opt_cut(rho, affine, tol);
      crho -= gen.d_data.d_prob * gen.crho();
    }
    rho += crho / (1 + ret.d_tau.back());
  } while (crho > tol);

  return ret;
}


void Stagewise::solve(int stage, int node, bool affine, bool lp, bool force)
{
  Master &master = d_masters[stage][node];
  master.solve_lp();

  if (lp)
    return;

  while (not master.integer())
  {
    Cut fenchel_cp = fenchel_cut(stage, node, affine);
    if (not add_cp(fenchel_cp, stage, node))    // mip could not be solved using cutting planes
    {
      if (force)
        master.solve_mip();                        // use Gurobi
      return;
    }

    master.solve_lp();
  }

  d_masters[stage][node].solve_lp();
}


Cut Stagewise::sddp_cut(int stage, Solution const &sol)
{
  Cut ret(nvars(stage));
  vmaster &subs = d_masters[stage + 1];

  for (int child = 0; child != subs.size(); ++child)
  {
    subs[child].update(sol);
    ret += d_stages[stage + 1][child].d_prob * subs[child].opt_cut();
  }
  return ret;
}

Cut Stagewise::fenchel_cut(int stage, int node, bool affine, double tol)
{
  return d_fenchel[stage][node].feas_cut(d_masters[stage][node].forward(), affine, tol);
}

void Stagewise::init_enums(int stage, const Solution &sol)
{
  for (int child = 0; child != outcomes(stage); ++child)
  {
    Master &sub = d_masters[stage][child];
    sub.update(sol);
    sub.solve_mip();

    Enumerator &gen = d_enumerators[stage][child];
    gen.set_mp(sol);
    gen.add_point(sub.forward(), false, true);
  }
}

bool Stagewise::add_cp(Cut &cut, int stage, int node, double tol)
{
  if (cut.d_tau.back() > 0)
    cut.scale();

  return d_masters[stage][node].add_cut(cut, tol);
}

void Stagewise::add_cut(Cut &cut, int stage)        // assumes shared outer approximations
{
  if (cut.d_tau.back() > 0)
    cut.scale();

  for (Master &master : d_masters[stage])
    master.add(cut);

  for (Enumerator &fenchel : d_fenchel[stage])
    fenchel.add_cut(cut);

  for (int lvl = stage; lvl != stage + 2; ++lvl)
  {
    for (Enumerator &gen : d_enumerators[lvl])
      gen.add_cut(cut);
  }
}

vector<vpath> Stagewise::enumerate_paths(vector<vpath> paths)
{
  int stage = paths[0].size();

  vector<vpath> ret;
  int outcomes = d_stages[stage].size();
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

  if (stage != d_stages.size() - 1)
    return enumerate_paths(ret);

  return ret;
}


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

  d_stages[stage - 1].push_back(data);
}

vector<int> Stagewise::nvars(int stage) const
{
  vector<int> ret(stage + 1);
  for (int lvl = 0; lvl != stage + 1; ++lvl)
    ret[lvl] = d_stages[lvl][0].nvars();
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

int Stagewise::outcomes(int stage) const
{
  return d_stages[stage].size();
}