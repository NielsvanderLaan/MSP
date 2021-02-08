#include "nodedata.h"

void NodeData::to_box(bool wide)
{
  if (d_stage == 1)
    return;

  arma::uvec inds = wide ? arma::uvec{} : d_fixed_constrs;

  d_Amat = d_Amat.cols(inds);
  d_Bmat = d_Bmat.cols(inds);
  d_rhs = d_rhs.elem(inds);
  d_senses = d_senses.elem(inds);

  d_fixed_constrs = arma::linspace<arma::uvec>(0, ncons() - 1, ncons());
}