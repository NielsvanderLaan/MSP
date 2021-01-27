#include "master.h"

Cut Master::opt_cut()
{
  solve_lp();
  arma::vec pi = multipliers();

  size_t nlevels = d_data.d_stage - 1;
  size_t ncuts = d_cuts.size();

  vvec beta(nlevels);
  vdouble tau(nlevels, 0.0);
  double alpha = obj();

  for (int level = 0; level != nlevels; ++level)
    beta[level] = arma::vec(d_state.d_x[level].n_elem, arma::fill::zeros);

  arma::vec lambda(pi.memptr(), d_data.ncons(), false);   // constraint multipliers
  beta.back() = d_data.d_Bmat * lambda;

  double *mu = pi.memptr() + d_data.ncons();          // cut multipliers

  for (int level = 0; level != nlevels; ++level)
  {
    for (size_t cut = 0; cut != ncuts; ++cut)
    {
      beta[level] += mu[cut] * d_cuts[cut].d_beta[level];
      tau[level] += mu[cut] * d_cuts[cut].d_tau[level];
    }
    alpha += dot(beta[level], d_state.d_x[level]) + tau[level] * d_state.d_theta[level];
  }
  Cut cut {alpha, beta, tau};

  return Cut {alpha, beta, tau};
}