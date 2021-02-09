#include "stagewise.h"

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

Cut Stagewise::sddp_cut(int stage, Solution const &sol)
{
  Cut ret(nvars(stage));

  for (int child = 0; child != outcomes(stage + 1); ++child)
  {
    Master &sub = get_master(stage + 1, child);
    sub.update(sol);
    ret += d_stages[stage + 1][child].d_prob * sub.opt_cut();
  }

  return ret;
}

Cut Stagewise::fenchel_cut(int stage, int node, bool affine, double tol)
{
  return get_fenchel(stage, node).feas_cut(solution(stage, node), affine, tol);
}



void Stagewise::init_enums(int stage, int node, const Solution &sol)
{
  int future = stage + 1;
  for (int child = 0; child != outcomes(future); ++child)
  {
    Master &sub = get_master(future, child);

    sub.update(sol);
    sub.solve_mip();

    Enumerator &gen = get_enums(stage, node)[child];
    gen.set_mp(sol);
    gen.add_point(sub.forward(), false, true);
  }
}