#include "dbenders.h"

void dBenders::init_gomory(GRBEnv &env)
{
  size_t nstages = d_data.nstages();
  d_gomory.reserve(nstages);

  for (int stage = 0; stage != nstages; ++stage)
  {
    vgom goms;
    size_t n_outcomes = d_data.outcomes(stage);
    goms.reserve(n_outcomes);

    bool leaf = stage == nstages - 1;
    for (int out = 0; out != n_outcomes; ++out)
      goms.emplace_back(Gomory{env,
                               node_data(stage, out),
                               leaf});

    d_gomory.emplace_back(move(goms));
  }
}