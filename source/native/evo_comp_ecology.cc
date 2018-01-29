// This is the main function for the NATIVE version of this project.

#include <iostream>
#include "../ecology_world.h"

#ifndef EMP_TRACK_MEM
#define EMP_TRACK_MEM
#endif

int main(int argc, char* argv[])
{
  
  BoxConfig config;
  auto args = emp::cl::ArgManager(argc, argv);
  if (args.ProcessConfigOptions(config, std::cout, "BoxConfig.cfg", "Box-macros.h") == false) exit(0);
  if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

  emp::Random rnd(config.SEED());

  if (config.PROBLEM() == "box") {
    using ORG_TYPE = emp::vector<double>;
    EcologyWorld<ORG_TYPE> world(rnd);
    world.Setup(config);
    world.Run();
  } else {
    using ORG_TYPE = emp::AvidaGP;
    EcologyWorld<ORG_TYPE> world(rnd);
    world.Setup(config);
    world.Run();
  }

}
