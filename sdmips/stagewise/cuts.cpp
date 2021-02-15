#include "stagewise.h"

Cut Stagewise::scaled_cut(int stage, int node, const Solution &sol, bool affine, double tol)
{
  Cut ret;
  vector<int> path_nvars = nvars(stage);

  double rho = sol.theta_n();
  double crho;

  //init_enums(stage, node, sol);

  v_enum enums = raw_enums(stage, node, sol);

  do
  {
    ret = Cut(path_nvars);
    crho = -rho;

    //for (Enumerator &gen : get_enums(stage, node))
    for (Enumerator &gen : enums)
    {
      ret += gen.d_data.d_prob * gen.opt_cut(rho, affine, tol);
      crho -= gen.d_data.d_prob * gen.crho();
    }
    rho += crho / (1 + ret.d_tau.back());
  } while (crho > tol);

  return ret;
}

Cut Stagewise::shared_scaled_cut(int stage, const Solution &sol, bool affine, double tol)
{
  assert(d_depth == 0);
  return scaled_cut(stage, 0, sol, affine, tol);
}

Cut Stagewise::sddp_cut(int stage, Solution const &sol)
{
  Cut ret(nvars(stage));

  int future = stage + 1;
  for (int child = 0; child != outcomes(future); ++child)
  {
    Master &sub = get_master(future, child);
    sub.update(sol);
    ret += d_stages[future][child].d_prob * sub.opt_cut();
  }

  return ret;
}

void Stagewise::init_enums(int stage, int node, Solution const &sol)
{
  v_enum &enumerators = get_enums(stage, node);

  vector<int> childs = children(stage, node);
  for (int idx = 0; idx != childs.size(); ++idx)
  {
    Master &sub = get_master(stage + 1, childs[idx]);

    sub.update(sol);
    sub.solve_mip();

    enumerators[idx].set_mp(sol);
    enumerators[idx].add_point(sub.forward(), true);
  }
}


v_enum Stagewise::raw_enums(int stage, int node, Solution const &sol)
{
  int future = stage + 1;

  v_enum enumerators;
  enumerators.reserve(outcomes(future));

  vpath path(stage + 2);
  iota(path.begin(), path.end(), 0);

  vector<NodeData> nodes;
  nodes.reserve(path.size());
  int box = max(stage - d_depth + 1, 0);
  for (int lvl = 0; lvl < box; ++lvl)
  {
    nodes.push_back(d_stages[lvl].front());
    nodes.back().to_box();
  }

  vector<int> sub_path = tail(stage, node);
  for (int depth = 0; depth != sub_path.size(); ++depth)
    nodes.push_back(d_stages[box + depth][sub_path[depth]]);

  nodes.resize(nodes.size() + 1);

  bool leaf = (future == d_stages.size() - 1);
  for (NodeData const &data : d_stages[future])
  {
    nodes.back() = data;
    enumerators.emplace_back(Enumerator{nodes, path, path.size() - 1, leaf, d_env});
  }

    // add cuts of parent
  Master &mp = get_master(stage, node);
  for (Enumerator &gen : enumerators)
  {
    for (Cut const &cut : mp.d_cuts)
      gen.add_cut(cut);
  }
    // add contemporary cuts
  vector<int> childs = children(stage, node);
  for (int idx = 0; idx != childs.size(); ++idx)
  {
    for (Cut const &cut : get_master(future, childs[idx]).d_cuts)
      enumerators[idx].add_cut(cut);

    Master &sub = get_master(future, childs[idx]);

    sub.update(sol);
    sub.solve_mip();

    enumerators[idx].set_mp(sol);
    enumerators[idx].add_point(sub.forward(), true);
  }

  return enumerators;
}