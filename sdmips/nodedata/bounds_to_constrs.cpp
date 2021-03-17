#include "nodedata.h"

void NodeData::bnds_to_cons()
{
  arma::mat Amat(d_Amat);
  arma::mat Bmat(d_Bmat);


  for (size_t var = 0; var != nvars(); ++var)
  {
    double lb = d_lb[var];
    double ub = d_ub[var];

    if (lb > 0)
    {
      Amat.insert_cols(Amat.n_cols, 1);
      Amat.col(Amat.n_cols - 1)[var] = 1.0;
      Bmat.insert_cols(Bmat.n_cols, 1);
      d_senses.resize(d_senses.n_elem + 1);
      d_senses.back() = GRB_GREATER_EQUAL;

      d_rhs.resize(d_rhs.n_elem + 1);
      d_rhs.back() = lb;

      int idx = Bmat.n_cols - 1;

      d_fixed_constrs.resize(d_fixed_constrs.n_elem + 1);
      d_fixed_constrs.back() = idx;

      d_box_constrs.resize(d_box_constrs.n_elem + 1);
      d_box_constrs.back() = idx;
    }

    if (ub < 1e20)
    {
      Amat.insert_cols(Amat.n_cols, 1);
      Amat.col(Amat.n_cols - 1)[var] = 1.0;
      Bmat.insert_cols(Bmat.n_cols, 1);
      d_senses.resize(d_senses.n_elem + 1);
      d_senses.back() = GRB_LESS_EQUAL;

      d_rhs.resize(d_rhs.n_elem + 1);
      d_rhs.back() = ub;

      int idx = Bmat.n_cols - 1;

      d_fixed_constrs.resize(d_fixed_constrs.n_elem + 1);
      d_fixed_constrs.back() = idx;

      d_box_constrs.resize(d_box_constrs.n_elem + 1);
      d_box_constrs.back() = idx;
    }

  }

  d_Amat = arma::sp_mat(Amat);
  d_Bmat = arma::sp_mat(Bmat);

  fill_n(d_lb.begin(), nvars(), 0.0);
  fill_n(d_ub.begin(), nvars(), GRB_INFINITY);
}