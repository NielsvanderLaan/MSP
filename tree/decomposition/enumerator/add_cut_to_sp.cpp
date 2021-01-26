#include "enumerator.h"

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