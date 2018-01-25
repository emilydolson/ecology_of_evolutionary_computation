// This is the main function for the NATIVE version of this project.

#include <iostream>
#include "../ecology_world.h"

#ifndef EMP_TRACK_MEM
#define EMP_TRACK_MEM
#endif

int main(int argc, char* argv[])
{
  EcologyWorld world;
  auto args = emp::cl::ArgManager(argc, argv);
  if (args.ProcessConfigOptions(world.config, std::cout, "BoxConfig.cfg", "Box-macros.h") == false) exit(0);
  if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

  world.Setup();
  world.Run();
}
