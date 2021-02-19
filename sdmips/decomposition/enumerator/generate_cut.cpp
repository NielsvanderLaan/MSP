#include "enumerator.h"

Cut Enumerator::fdecom(double tol, bool affine, bool reset)
{
  if (reset) clear();
  set_mp(false);              // relax mp
  set_sp(false);              // relax sp

  Cut cut;
  double diff = GRB_INFINITY;
  bool mp_strike = false;
  bool sp_strike = false;

  while (true)
  {
    try
    {
      cut = solve_mp(affine, 1e8);  // candidate
      set_sub(cut);
      solve_sp();

      if (mp_violation() >= max(diff - 1e-6, tol))
        throw mp_exception{};

      diff = cut.d_alpha - sub_val();
      if (diff < tol)
        break;
      add_point(point());

      if (mp_strike) set_mp(false);
      mp_strike = false;

      if (sp_strike) set_sp(false);
      sp_strike = false;
    } catch (mp_exception &e)
    {
      if (not mp_strike)
      {
        mp_strike = true;
        set_mp(true);
        continue;
      }
      if (not reset)
        return fdecom(tol, affine, true);
      cout << "Fdecom: unrecoverable numerical difficulties, gap: " << cut.d_alpha - sub_bound() << '\n';
      break;
    } catch (sp_exception &e)
    {
      if (not sp_strike)
      {
        sp_strike = true;
        set_sp(true);
        continue;
      }

      cout << "sp_status: " << sp_status() << '\n';
      cut.d_alpha = -GRB_INFINITY;
      return cut;
    }
  }

  cut.d_alpha = sub_bound();
  assert(cut.depth() > 0);

  return cut;
}

Cut Enumerator::feas_cut(Solution const &sol, bool affine, double tol)
{
  set_mp(sol);
  set_bounds(affine, 1e2);
  Cut ret = fdecom(tol, affine);

  ret.d_tau.back() -= 1;
  ret.d_feas = true;
  return ret;
}

Cut Enumerator::opt_cut(double rho, bool affine, double tol)       // bool affine
{
  set_tau_bounds(affine, GRB_INFINITY);
  set_rho(rho);
  return fdecom(tol, affine);
}

