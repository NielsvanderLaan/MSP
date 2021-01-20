#include "structs.h"

double compute_lhs(Cut &cut, Solution &sol)
{
  int range = sol.depth() - cut.depth();
  assert(range >= 0);

  return (1 + cut.d_tau[0]) * sol.d_theta[range] + arma::dot(cut.d_beta[0], sol.d_x[range]);
}

double compute_rhs(Cut &cut, Solution &sol)
{
  int range = sol.depth() - cut.depth();
  assert(range >= 0);

  double ret = cut.d_alpha;

  for (size_t idx = 1; idx < cut.depth(); ++idx)
  {
    ret -= arma::dot(cut.d_beta[idx], sol.d_x[idx + range]);
    ret -= cut.d_tau[idx] * sol.d_theta[idx + range];
  }

  return ret;
}

double scaled_rhs(Cut &cut, Solution &sol)
{
  int range = sol.depth() - cut.depth();
  assert(range >= 0);

  double numerator = cut.d_alpha;

  for (size_t idx = 0; idx < cut.depth(); ++idx)
    numerator -= arma::dot(cut.d_beta[idx], sol.d_x[idx + range]);

  for (size_t idx = 1; idx < cut.depth(); ++idx)
    numerator -= cut.d_tau[idx] * sol.d_theta[idx + range];

  return numerator / (1 + cut.d_tau[0]);
}
















