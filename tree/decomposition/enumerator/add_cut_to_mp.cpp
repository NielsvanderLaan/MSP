#include "enumerator.h"

void Enumerator::add_cut_to_mp(Cut &cut)
{
  int range = d_theta.size() - cut.depth();

  GRBConstr *mp_cons = d_mp.getConstrs();

  for (size_t con = 0; con != d_points.size(); ++con)
  {
    Solution &point = d_points[con];
    double diff = max(0.0, scaled_rhs(cut, point) - point.d_theta[range]);
    point.d_theta[range] += diff;
    GRBConstr &mp_cut = mp_cons[con];

    if (range > 0)
      d_mp.chgCoeff(mp_cut, d_tau[range - 1], point.d_theta[range - 1]);
    else
      mp_cut.set(GRB_DoubleAttr_RHS, mp_cut.get(GRB_DoubleAttr_RHS) + diff);
  }

  delete[] mp_cons;
  d_mp.update();
}