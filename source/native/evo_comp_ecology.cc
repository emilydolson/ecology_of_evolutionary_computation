// This is the main function for the NATIVE version of this project.

#include <iostream>
#include "../ecology_world.h"

#ifndef EMP_TRACK_MEM
#define EMP_TRACK_MEM
#endif

int main()
{
  EcologyWorld world("testcases/count-odds.csv");
  world.Setup();
  world.Run();
}
