#include "tree.h"

void Tree::SND(bool affine)
{
  while (true)
  {
    forward(affine, false);
    cout << d_masters[0].mip_obj() << '\n';   // todo: lp or mip

    if (not backward(affine))    // no improvement possible
      break;
  }
}