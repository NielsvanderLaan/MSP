#include "benders.h"

Cut Benders::compute_cut(Family type, Solution const &sol, int stage, int node)
{
  switch (type)
  {
    case SDDP:
      return sddp_cut(stage, sol);
    case LR:
      return scaled_cut(stage, node, sol, true);
    case SC:
      return scaled_cut(stage, node, sol, false);
    case LBDA_ZEROS:
      return lbda_cut(stage, sol, ZEROS);
    case LBDA_RC:
      return lbda_cut(stage, sol, RECURSIVE);

    default:
      cout << "unknown cut type\n";
      exit(1);
  }
}

Cut Benders::scaled_cut(int stage, int node, const Solution &sol, bool affine, double tol)
{
  Cut ret;
  vector<int> path_nvars = d_data.nvars(stage);

  double rho = sol.theta_n();
  double crho;

  v_enum &gens = init_enums(stage, node, sol);

  vector<double> probs = d_data.probs( stage + 1);

  do
  {
    ret = Cut(path_nvars);
    crho = -rho;

    for (size_t child = 0; child != gens.size(); ++child)
    {
      ret += probs[child] * gens[child].opt_cut(rho, affine, tol);
      crho -= probs[child] * gens[child].crho();
    }
    if (affine) break;
    rho += crho / (1 + ret.d_tau.back());
  } while (crho > tol);

  return ret;
}

Cut Benders::sddp_cut(int stage, Solution const &sol)
{
  Cut ret(d_data.nvars(stage));

  int future = stage + 1;
  vector<double> probs = d_data.probs(future);
  for (int child = 0; child != d_data.outcomes(future); ++child)
  {
    Master &sub = get_master(future, child);
    sub.update(sol);
    ret += probs[child] * sub.opt_cut();
  }

  return ret;
}

Cut Benders::lbda_cut(int stage, const Solution &sol, Alpha type, arma::vec alpha)
{
  int future = stage + 1;

  if (type == ZEROS)
    alpha = arma::vec(node_data(future, 0).ncons(), arma::fill::zeros);

  Cut ret(d_data.nvars(stage));
  vector<double> probs = d_data.probs(future);

  for (int out = 0; out != d_data.outcomes(future); ++out)
  {
    Master &sub = get_master(future, out);
    Gomory &gom = get_gomory(future, out);

    sub.update(sol);
    sub.solve_lp();

    if (type == RECURSIVE)
      alpha = node_data(future, out).d_Bmat.t() * sol.d_x.back();

    gom.update(node_data(future, out).d_rhs - alpha, sub.vbasis(), sub.cbasis());
    /*
     * correctly takes into account deterministic constraints, but upper/lower bounds on
     * variables should be handled carefully (including theta: d_L).
     */
    //gom.update(node_data(future, out).d_rhs - alpha, sub.basis());

    gom.solve();

    auto lambda = sub.multipliers(false);
    ret.d_beta.back() += probs[out] * node_data(future, out).d_Bmat * lambda;
    ret.d_alpha += probs[out] * (gom.obj() + dot(lambda, alpha));
  }

  return ret;
}

v_enum &Benders::init_enums(int stage, int node, Solution const &sol)
{
  v_enum &gens = get_enums(stage, node);

  vector<int> childs = children(stage, node);
  for (int idx = 0; idx != childs.size(); ++idx)
  {
    Master &sub = get_master(stage + 1, childs[idx]);

    sub.update(sol);
    sub.solve_mip();

    gens[idx].set_mp(sol);
    gens[idx].add_point(sub.forward(), true);
  }

  return gens;
}
