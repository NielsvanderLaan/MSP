#ifndef MSP_FAMILIES_H
#define MSP_FAMILIES_H

#include<iostream>

using namespace std;

enum Family
{
    SDDP,
    LR,
    SC,
    DEFAULT
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
    default:
      return "unknown cut type\n";
  }
}

#endif //MSP_FAMILIES_H
