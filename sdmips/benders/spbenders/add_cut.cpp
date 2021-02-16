#include "spbenders.h"

void spBenders::add_cut(Cut &cut, int stage, vector<int> const &path)
{
  if (stage == 0)
    d_root.push_cut(cut);

  d_cuts[stage][master_idx(stage, path)].push_back(cut);
}

void spBenders::add_shared_cut(Cut &cut, int stage)
{
  if (stage == 0)
    d_root.push_cut(cut);

  for (outer_apx &apx : d_cuts[stage])
    apx.push_back(cut);
}