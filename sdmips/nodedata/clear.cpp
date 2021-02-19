#include "nodedata.h"

void NodeData::clear()
{
  d_Amat = arma::sp_mat(nvars(), 0);
  d_Bmat.clear();
  d_rhs.clear();
  d_senses.clear();



  d_fixed_constrs.clear();
  d_box_constrs.clear();
}