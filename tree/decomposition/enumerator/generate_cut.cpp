#include "enumerator.h"

Cut Enumerator::generate_cut(double rho, double tol)
{
  set_rho(rho);

  Cut cut;
  while (true)
  {
    solve_mp();
    cut = candidate();
    set_sub(cut);
    solve_sp();

    if (sub_val() >= cut.d_alpha - tol)
      break;

    add_point(point());
  }

  cut.d_alpha = sub_bound();
  return cut;
}