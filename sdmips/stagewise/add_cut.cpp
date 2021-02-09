#include "stagewise.h"

bool Stagewise::add_cp(Cut &cut, int stage, int node, double tol)
{
  if (cut.d_tau.back() > 0)
    cut.scale();

  return get_master(stage, node).add_cut(cut, tol);
}

void Stagewise::add_cut(Cut &cut, int stage, vector<int> const &path)
{
  if (cut.d_tau.back() > 0)
    cut.scale();

  int node = master_idx(stage, path);
  get_master(stage, node).add(cut);        // no sharing
  get_fenchel(stage, node).add_cut(cut);   // idem

  for (Enumerator &gen : get_enums(stage, node))
    gen.add_cut(cut);

  if (stage == 0)
    return;
  for (size_t parent : parents(stage, path))      // if depth = 0, then enum objects are shared and this is stupid{}
    get_enums(stage - 1, parent)[path.back()].add_cut(cut);

}
