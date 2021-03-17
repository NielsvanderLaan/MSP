#include "master.h"

void Master::solve_lp()
{
  d_lp->optimize();
  assert(d_lp->get(GRB_IntAttr_Status) == 2);
  d_x_n = lp_xvals();
  d_theta_n = lp_theta();
}

void Master::solve_mip()
{
  d_mip->optimize();
  assert(d_mip->get(GRB_IntAttr_Status) == 2);

  d_x_n = mip_xvals();
  d_theta_n = mip_theta();
}

arma::vec Master::mip_xvals()
{
  double *xvals = d_mip->get(GRB_DoubleAttr_X, d_mip_xvars.data(), d_mip_xvars.size());
  arma::vec ret(xvals, d_mip_xvars.size());
  delete[] xvals;
  return ret;
}

double Master::mip_theta()
{
  return d_mip_theta.get(GRB_DoubleAttr_X);
}

arma::vec Master::lp_xvals()
{
  double *xvals = d_lp->get(GRB_DoubleAttr_X, d_lp_xvars.data(), d_lp_xvars.size());
  arma::vec ret(xvals, d_lp_xvars.size());
  delete[] xvals;
  return ret;
}

double Master::lp_theta()
{
  return d_lp_theta.get(GRB_DoubleAttr_X);
}

Solution Master::forward()
{
  Solution sol = d_state;
  sol.extend(d_x_n, d_theta_n);
  return sol;
}

double Master::theta_n() const
{
  return d_theta_n;
}

bool Master::integer()
{
  return is_integer(lp_xvals(), d_data.d_types);
}

arma::vec Master::multipliers(bool cuts)
{
  int ncons = d_data.ncons();
  if (cuts)
    ncons += d_cuts.size();

  GRBConstr *cons = d_lp->getConstrs();

  double *pi = d_lp->get(GRB_DoubleAttr_Pi, cons, ncons);
  arma::vec ret(pi, ncons);

  delete[] cons;
  delete[] pi;
  return ret;
}

vector<int> Master::basis() const
{
  GRBmodel *cm = d_lp->Cmodel;
  int len = d_data.ncons() + d_cuts.size();

  int bhead[len];
  fill_n(bhead, len, -1);
  GRBgetBasisHead(cm, bhead);

  int Brow_inds[len];
  double Brow_vals[len];
  GRBsvec Brow{len, Brow_inds, Brow_vals};

  int inds[1];
  double vals[1] = {1.0};
  GRBsvec e_i{1, inds, vals};

  vector<int> ret;
  ret.reserve(len);

  for (int row = 0; row != len; ++row)
  {
    e_i.ind[0] = row;
    GRBBSolve(cm, &e_i, &Brow);

    for (int nz = 0; nz != Brow.len; ++nz)
    {
      if (Brow.ind[nz] < d_data.d_rn_constrs)
      {
        ret.push_back(bhead[row]);
        continue;
      }
    }
  }

  return ret;
}

vector<int> Master::vbasis() const
{
  int nvars = d_data.nvars();
  int *vb = d_lp->get(GRB_IntAttr_VBasis, d_lp_xvars.data(), nvars);

  vector<int> ret(nvars + 1);
  copy_n(vb, nvars, ret.begin());
  ret.back() = d_lp_theta.get(GRB_IntAttr_VBasis);

  return ret;
}

vector<int> Master::cbasis() const
{
  GRBConstr *cons = d_lp->getConstrs();
  int len = d_data.ncons() + d_cuts.size();
  int *cb = d_lp->get(GRB_IntAttr_CBasis, cons, len);

  vector<int> ret(cb, cb + len);
  delete[] cb;
  delete[] cons;

  return ret;
}


double Master::lp_obj() const
{
  return d_lp->get(GRB_DoubleAttr_ObjBound);
}


double Master::mip_obj() const
{
  return d_mip->get(GRB_DoubleAttr_ObjBound);
}

void Master::set_rho(double rho)
{
  double diff = rho - d_state.d_theta.back() ;
  d_state.d_theta.back() = rho;
  GRBConstr *cons = d_lp->getConstrs();
  double *rhs = d_lp->get(GRB_DoubleAttr_RHS, cons + d_data.ncons(), d_cuts.size());

  for (size_t cut = 0; cut != d_cuts.size(); ++cut)
  {
    vdouble const &tau = d_cuts[cut].d_tau;
    rhs[cut] -= diff * tau[tau.size() - 2];
  }

  d_lp->set(GRB_DoubleAttr_RHS, cons + d_data.ncons(), rhs, d_cuts.size());
  d_lp->update();

  delete[] cons;
  delete[] rhs;
}