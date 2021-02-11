#include "stagewise.h"

#include <fstream>

Cut Stagewise::scaled_cut(int stage, int node, const Solution &sol, bool affine, double tol)   // TODO: sharing
{
  Cut ret;
  vector<int> path_nvars = nvars(stage);

  double rho = sol.theta_n();
  double crho;

  init_enums(stage, node, sol);

  do
  {
    ret = Cut(path_nvars);
    crho = -rho;

    for (Enumerator &gen : get_enums(stage, node))
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

Cut Stagewise::fenchel_cut(int stage, int node, bool affine, double tol)
{
  return get_fenchel(stage, node).feas_cut(solution(stage, node), affine, tol);
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