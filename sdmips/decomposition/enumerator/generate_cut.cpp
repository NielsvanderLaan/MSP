#include "enumerator.h"

Cut Enumerator::fdecom(double tol, bool reset)
{
  d_mp->write("mp_before.lp");
  if (reset) clear();

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
      cout << "adding direction\n";
      add_point(direction(), true);
      continue;
    }

    diff = cut.d_alpha - sub_val();
    if (diff < tol)
      break;

    if (mp_violation() > max(diff - 1e-6, tol))
      goto numerical;

    add_point(point());

    first_strike = false;
    set_mp(false);
    continue;

numerical:                // numerical difficulties occured
    cout << "mp status: " << mp_status() << ". violation: " << mp_violation() << '\n';
    cout << d_data.d_stage << '\n';
    d_mp->write("mp.lp");
    exit(8);
    if (not first_strike)
    {
      first_strike = true;
      set_mp(true);
      continue;
    }
    if (not reset)
      return fdecom(tol, true);

    cerr << "Fdecom: unrecoverable numerical difficulties.\n";
    break;
  }

  set_mp(false);              // relax mp
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

