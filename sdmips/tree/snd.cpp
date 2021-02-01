#include "tree.h"

void Tree::SND()
{
  // TODO: time limit, print more info

  while (true)
  {
    forward(false);
    cout << d_masters[0].obj() << '\n';

    if (not backward())    // no improvement possible
      break;
  }
}