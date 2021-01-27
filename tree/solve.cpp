#include "tree.h"

void Tree::solve()
{
  // TODO: time limit, print more info

  while (true)
  {
    forward(true);
    cout << d_masters[0].obj() << '\n';

    if (not backward())    // no improvement possible
      break;
  }
}