#include "spbenders.h"

void spBenders::add_cut(Cut &cut, int stage, int node)
{
  if (stage == 0)
    d_root.push_cut(cut);

  d_nodal_apx[stage][node].push_back(cut);
}

void spBenders::add_shared_cut(Cut &cut, int stage)
{
  if (stage == 0)
    d_root.push_cut(cut);

  d_stage_apx[stage].push_back(cut);
}