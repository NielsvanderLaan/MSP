#include "enumerator.h"

void Enumerator::add_cut_to_mp(Cut &cut)
{
  int depth = cut.depth();

  GRBConstr *mp_cons = d_mp.getConstrs();

  for (size_t con = 0; con != d_points.size(); ++con)
  {
    Solution &point = d_points[con];
    double diff = max(0.0, scaled_rhs(cut, point) - point.d_theta[depth - 1]);
    point.d_theta[depth - 1] += diff;
    GRBConstr &mp_cut = mp_cons[con];

    if (depth < d_data.d_stage)
      d_mp.chgCoeff(mp_cut, d_tau[depth - 1], -point.d_theta[depth - 1]);
    else
      mp_cut.set(GRB_DoubleAttr_RHS, mp_cut.get(GRB_DoubleAttr_RHS) + diff);
  }

  delete[] mp_cons;
  d_mp.update();
}