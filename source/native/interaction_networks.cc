#include "interaction_networks.h"
#include "tools/string_utils.h"
#include "config/ArgManager.h"
#include "config/command_line.h"

EMP_BUILD_CONFIG( InteractionConfig,
  VALUE(REPS_PER_CONDITION, size_t, 10, "How many replicates per condition?"),
  VALUE(POP_SIZE_START, size_t, 10, "Lowest population size"),
  VALUE(POP_SIZE_END, size_t, 100, "Highest population size"),
  VALUE(POP_SIZE_INC, size_t, 10, "Increment of population sizes"),
  VALUE(N_TRAITS_START, size_t, 5, "Lowest number of traits"),
  VALUE(N_TRAITS_END, size_t, 50, "Highest number of traits"),
  VALUE(N_TRAITS_INC, size_t, 5, "Increment of number of traits"),
  VALUE(SIGMA_SHARE_START, double,0, "Lowest sigma share"),
  VALUE(SIGMA_SHARE_END, double, 10, "Highest sigma share"),
  VALUE(SIGMA_SHARE_INC, double, 100, "Increment of sigma share"),
  VALUE(ALPHA_START, double,0, "Lowest alpha"),
  VALUE(ALPHA_END, double, 2, "Highest alpha"),
  VALUE(ALPHA_INC, double, .2, "Increment of alpha"),
  VALUE(COST_START, double,0, "Lowest COST"),
  VALUE(COST_END, double, 10, "Highest COST"),
  VALUE(COST_INC, double, 1, "Increment of COST"),
  VALUE(CF_START, double, 0.1, "Lowest CF"),
  VALUE(CF_END, double, 0.1, "Highest CF"),
  VALUE(CF_INC, double, 1, "Increment of CF"),
  VALUE(NICH_WIDTH_START, double, 0, "Lowest NICH_WIDTH"),
  VALUE(NICH_WIDTH_END, double, 10, "Highest NICH_WIDTH"),
  VALUE(NICH_WIDTH_INC, double, 1, "Increment of NICH_WIDTH"),
  VALUE(MAX_SCORE_START, double, 1, "Lowest MAX_SCORE"),
  VALUE(MAX_SCORE_END, double, 101, "Highest MAX_SCORE"),
  VALUE(MAX_SCORE_INC, double, 10, "Increment of MAX_SCORE")

)


int main(int argc, char* argv[]) {
    InteractionConfig config;
    config.Read("Interaction.cfg");
    auto args = emp::cl::ArgManager(argc, argv);
    if (args.ProcessConfigOptions(config, std::cout, "Interaction.cfg", "Interaction-macros.h") == false) exit(0);
    if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

    Controller c;    

    for (size_t n_traits = config.N_TRAITS_START(); n_traits <= config.N_TRAITS_END(); n_traits += config.N_TRAITS_INC()) {
        c.SetNTraits(n_traits);
        std::cout << "n traits" << n_traits << std::endl;
        for (size_t pop_size = config.POP_SIZE_START(); pop_size <= config.POP_SIZE_END(); pop_size += config.POP_SIZE_INC()) {
            c.SetPopSize(pop_size);
            std::cout << pop_size << std::endl;
            for (size_t rep = 0; rep < config.REPS_PER_CONDITION(); rep++) {
                std::cout << "rep:" << rep << std::endl;
                c.Regenerate();
            }
        }
    }
}