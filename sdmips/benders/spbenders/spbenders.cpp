#include "spbenders.h"

spBenders::spBenders(GRBEnv &env, Stagewise &data, int depth)
        :
        Benders(env, data, depth),
        d_stage_apx(d_data.nstages() - 1),
        d_nodal_apx(d_data.nstages() - 1),
        d_root(node_data(0, 0), false, env)
{
  for (int stage = 0; stage != d_nodal_apx.size(); ++stage)
  {
    size_t n_nodes = 1;

    for (size_t lvl = 0; lvl != min(d_depth, stage + 1); ++lvl)
      n_nodes *= d_data.outcomes(stage - lvl);

    d_nodal_apx[stage] = stage_apx(n_nodes);
  }
}