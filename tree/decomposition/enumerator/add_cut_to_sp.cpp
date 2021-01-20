#include "enumerator.h"

void Enumerator::add_cut_to_sp(Cut &cut)
{
  int depth = cut.depth();
  int range = d_theta.size() - depth;

  GRBLinExpr lhs = d_theta[range];
  for (size_t lvl = 0; lvl != depth; ++lvl)
  {
    lhs += cut.d_tau[lvl] * d_theta[range + lvl];
    vvar &x_n = d_x[range + lvl];
    lhs.addTerms(cut.d_beta[lvl].memptr(),x_n.data(),  x_n.size());
  }
  d_sp.addConstr(lhs >= cut.d_alpha);

  d_sp.update();
}