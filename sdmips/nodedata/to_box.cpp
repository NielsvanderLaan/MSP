#include "nodedata.h"

NodeData NodeData::to_box() const
{
  /*
  NodeData copy = *this;
  if (d_stage == 1)
    return copy;

  copy.d_Amat = d_Amat.cols(d_fixed_constrs);
  copy.d_Bmat = d_Bmat.cols(d_fixed_constrs);
  copy.d_rhs = d_rhs.elem(d_fixed_constrs);
  copy.d_senses = d_senses.elem(d_fixed_constrs);

  copy.d_fixed_constrs = arma::linspace<arma::uvec>(0,
                                                    copy.ncons() - 1,
                                                    copy.ncons());
  return copy;
   */
  return *this;
}