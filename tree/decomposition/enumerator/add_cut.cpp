#include "enumerator.h"

void Enumerator::add_cut(Cut &cut)
{
  add_cut_to_sp(cut);
  add_cut_to_mp(cut);
}
