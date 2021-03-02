#ifndef MSP_FAMILIES_H
#define MSP_FAMILIES_H

#include<iostream>

using namespace std;

enum Family
{
  SDDP,
  LR,
  SC,
  LBDA_ZEROS,
  LBDA_RC,
  DEFAULT
};

enum Alpha
{
  ZEROS,
  RECURSIVE
  // USER
  // RANDOM
};

static string to_string(Family type)
{
  switch (type)
  {
    case SDDP:
      return "SDDP cuts";
    case LR:
      return "Lagrangian cuts";
    case SC:
      return "Scaled cuts";
    case LBDA_ZEROS:
      return "LBDA (alpha = 0) cuts";
    case LBDA_RC:
      return "LBDA (recursive) cuts";
    default:
      return "unknown cut type\n";
  }
}

#endif //MSP_FAMILIES_H
