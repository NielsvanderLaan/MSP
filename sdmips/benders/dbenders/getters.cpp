#include "dbenders.h"

Master &dBenders::get_master(int stage, int node)
{
  return d_nodes[stage][node].first;
}

v_enum &dBenders::get_enums(int stage, int node)
{
  return *d_nodes[stage][node].second;
}

vector<outer_apx> dBenders::export_cuts()
{
  assert(d_depth == 0);
  vector<outer_apx> ret(d_data.nstages() - 1);

  for (int stage = 0; stage != ret.size(); ++stage)
  {
    Master &mp = get_master(stage, 0);
    ret[stage].reserve(mp.d_cuts.size());
    for (Cut const &cut : mp.d_cuts)
      ret[stage].push_back(cut);
  }

  return ret;
}
