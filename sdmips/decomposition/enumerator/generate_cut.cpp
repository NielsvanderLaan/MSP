#include "enumerator.h"

Cut Enumerator::fdecom(double tol, bool reset)
{
  if (reset) clear();
  set_mp(false);              // relax mp

  Cut cut;
  bool first_strike = false;

  while (true)
  {
    solve_mp();
    double diff;
    if (mp_status() != 2)
      goto numerical;

    cut = candidate();
    set_sub(cut);
    solve_sp();

    if (sp_status() != 2)
    {
      add_point(direction(), true);
      continue;
    }

    diff = cut.d_alpha - sub_val();
    if (diff < tol)
      break;

    if (mp_violation() >= max(diff - 1e-6, tol))
      goto numerical;

    add_point(point(), false);

    if (first_strike)
      set_mp(false);
    first_strike = false;

    continue;

numerical:                // numerical difficulties occured
    //cout << "mp status: " << mp_status() << ". violation: " << mp_violation() << '\n';
    if (not first_strike)
    {
      first_strike = true;
      set_mp(true);
      continue;
    }
    if (not reset)
      return fdecom(tol, true);

    cerr << "Fdecom: unrecoverable numerical difficulties, gap: " << cut.d_alpha - sub_bound() << '\n';
    break;
  }

  cut.d_alpha = sub_bound();
  assert(cut.depth() > 0);

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

