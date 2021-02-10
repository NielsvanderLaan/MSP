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

  if (stage == 0)     // root problem --> no parents
    return;

  for (size_t parent : parents(stage, path))
    get_enums(stage - 1, parent)[path[stage]].add_cut(cut);
}

void Stagewise::add_shared_cut(Cut &cut, int stage)
{
  if (cut.d_tau.back() > 0)
    cut.scale();

  for (size_t node = 0; node != d_nodes[stage].size(); ++node)
  {
    get_master(stage, node).add(cut);
    get_fenchel(stage, node).add_cut(cut);
  }

  for (Enumerator &gen : get_enums(stage, 0))             // enumerators are shared
    gen.add_cut(cut);

  if (stage == 0)
    return;

  for (Enumerator &gen : get_enums(stage - 1, 0))   // contemporary cut
    gen.add_cut(cut);
}