#include "benders.h"

void Benders::import_cuts(vector<outer_apx> stage_apx)
{
  int stage = 0;
  for (outer_apx &apx : stage_apx)
  {
    for (Cut &cut : apx)
      add_shared_cut(cut, stage);

    ++stage;
  }
}