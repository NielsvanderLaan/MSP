#include "dbenders.h"

void dBenders::add_cut(Cut &cut, int stage, int node)
{
  get_master(stage, node).push_cut(cut);

  for (Enumerator &gen : get_enums(stage, node))
    gen.add_cut(cut);

  if (stage == 0)     // root problem --> no parents
    return;

  int out = outcome(stage, node);
  for (size_t parent : parents(stage, node))
    get_enums(stage - 1, parent)[out].add_cut(cut);

}

void dBenders::add_shared_cut(Cut &cut, int stage)
{
  for (size_t node = 0; node != d_nodes[stage].size(); ++node)
    get_master(stage, node).push_cut(cut);


  for (Enumerator &gen : get_enums(stage, 0))             // enumerators are shared
    gen.add_cut(cut);

  if (stage == 0)
    return;

  for (Enumerator &gen : get_enums(stage - 1, 0))   // contemporary cut
    gen.add_cut(cut);
}