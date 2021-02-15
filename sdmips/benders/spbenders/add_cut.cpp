#include "spbenders.h"

void spBenders::add_cut(Cut &cut, int stage, vector<int> const &path)
{
  d_cuts[stage][master_idx(stage, path)].push_back(cut);
}

void spBenders::add_shared_cut(Cut &cut, int stage)
{
  for (outer_apx &apx : d_cuts[stage])
    apx.push_back(cut);
}