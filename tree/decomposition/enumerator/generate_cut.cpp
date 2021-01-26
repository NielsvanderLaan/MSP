#include "enumerator.h"

Cut Enumerator::fdecom(double tol)
{
  Cut cut;
  while (true)
  {
    solve_mp();
    cut = candidate();
    set_sub(cut);
    solve_sp();

    if (sub_val() > cut.d_alpha - tol)
      break;

    add_point(point());
  }

  cut.d_alpha = sub_bound();
  return cut;
}

Cut Enumerator::feas_cut(Solution const &sol, double tol)
{
  set_mp(sol);
  Cut ret = fdecom(tol);

  ret.d_tau.back() -= 1;
  ret.d_feas = true;
  return ret;
}

Cut Enumerator::opt_cut(double rho, double tol)
{
  set_rho(rho);
  return fdecom(tol);
}

