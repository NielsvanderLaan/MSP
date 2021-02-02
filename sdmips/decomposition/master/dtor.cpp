#include "master.h"

Master::~Master()
{
  if (d_mip)
    delete d_mip;

  if (d_lp)
    delete d_lp;
}