#include "enumerator.h"

void Enumerator::add_cut_to_mp(Cut &cut)
{
  int idx = cut.depth() - 1;

  GRBConstr *mp_cons = d_mp.getConstrs();

  for (size_t con = 0; con != d_points.size(); ++con)
  {
    Solution &point = d_points[con];
    double diff = max(0.0, scaled_rhs(cut, point) - point.d_theta[idx]);
    point.d_theta[idx] += diff;

    GRBConstr &mp_cut = mp_cons[con];
    if (idx < d_tau.size())
      d_mp.chgCoeff(mp_cut, d_tau[idx], -point.d_theta[idx]);
    else
      mp_cut.set(GRB_DoubleAttr_RHS, mp_cut.get(GRB_DoubleAttr_RHS) + diff);
  }

  delete[] mp_cons;
  d_mp.update();
}