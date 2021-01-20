#include "master.h"

Cut Master::compute_cut(Solution &sol)
{
  update(sol);
  optimize();
  arma::vec pi = multipliers();

  size_t nlevels = sol.depth();
  size_t ncuts = d_cuts.size();

  vvec beta(nlevels);
  vdouble tau(nlevels, 0.0);
  double alpha = obj();

  assert(nlevels == d_data.d_stage - 1);

  for (int level = 0; level != nlevels; ++level)
    beta[level] = arma::vec(sol.d_x[level].n_elem, arma::fill::zeros);

  arma::vec lambda(pi.memptr(), d_data.ncons(), false);   // constraint multipliers
  beta[0] = d_data.d_Bmat * lambda;

  double *mu = pi.memptr() + d_data.ncons();          // cut multipliers

  for (int level = 0; level != nlevels; ++level)
  {
    for (size_t cut = 0; cut != ncuts; ++cut)
    {
      beta[level] += mu[cut] * d_cuts[cut].d_beta[level];
      tau[level] += mu[cut] * d_cuts[cut].d_tau[level];
    }
    alpha += dot(beta[level], sol.d_x[level]) + tau[level] * sol.d_theta[level];
  }

  return Cut {alpha, beta, tau};
}