#include "enumerator.h"

void Enumerator::add_cut_to_sp(Cut &cut)
{
  int depth = cut.depth();

  GRBLinExpr lhs = d_theta[depth - 1];
  for (size_t stage = 0; stage != depth; ++stage)
  {
    lhs += cut.d_tau[stage] * d_theta[stage];
    vvar &x_n = d_x[stage];
    lhs.addTerms(cut.d_beta[stage].memptr(), x_n.data(),  x_n.size());
  }
  d_sp.addConstr(lhs >= cut.d_alpha);

  d_sp.update();
}