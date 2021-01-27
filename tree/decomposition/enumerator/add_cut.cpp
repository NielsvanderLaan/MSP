#include "enumerator.h"

void Enumerator::add_cut(Cut &cut)
{
  if (cut.d_feas)
    return;

  add_cut_to_sp(cut);
  add_cut_to_mp(cut);
}

void Enumerator::add_cut_to_sp(Cut &cut)
{
  int depth = cut.depth();

  GRBLinExpr lhs = d_theta[depth - 1];
  lhs.addTerms(cut.d_tau.data(), d_theta.data(), depth);

  for (size_t stage = 0; stage != depth; ++stage)
    lhs.addTerms(cut.d_beta[stage].memptr(), d_x[stage].data(), d_x[stage].size());

  d_sp.addConstr(lhs >= cut.d_alpha);

  d_sp.update();
}

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