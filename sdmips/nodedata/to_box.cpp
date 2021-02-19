#include "nodedata.h"

void NodeData::to_box(bool wide)
{
  if (d_stage == 1)
    return;

  arma::uvec inds = wide ? d_box_constrs : d_fixed_constrs;

  d_Amat = d_Amat.cols(inds);
  d_Bmat = d_Bmat.cols(inds);
  d_rhs = d_rhs.elem(inds);
  d_senses = d_senses.elem(inds);

  vector<size_t> box_constrs;
  for (int con : d_box_constrs)
    box_constrs.push_back(distance(d_fixed_constrs.begin(),
                                   find(d_fixed_constrs.begin(), d_fixed_constrs.end(), con)));
  d_box_constrs = arma::uvec(box_constrs);

  d_fixed_constrs = arma::linspace<arma::uvec>(0, ncons() - 1, ncons());
}