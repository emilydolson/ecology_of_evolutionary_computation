#include "config/ArgManager.h"
#include "base/vector.h"
#include "tools/Random.h"
#include "tools/math.h"
#include "tools/set_utils.h"
#include "Evo/World.h"
#include "Evo/Resource.h"
#include "hardware/AvidaGP.h"
#include "TestcaseSet.h"
#include "data/DataNode.h"
#include "control/Signal.h"
#include "tools/map_utils.h"
#include <map>

EMP_BUILD_CONFIG( BoxConfig,
  GROUP(DEFAULT, "Default settings for box experiment"),
  VALUE(SEED, int, 0, "Random number seed (0 for based on time)"),
  VALUE(POP_SIZE, uint32_t, 1000, "Number of organisms in the popoulation."),
  VALUE(START_POP_SIZE, uint32_t, 1, "Number of organisms to initialize population with."),
  VALUE(UPDATES, uint32_t, 1000, "How many generations should we process?"),
  VALUE(SELECTION, std::string, "LEXICASE", "What selection scheme should we use?"),
//   VALUE(N_NEUTRAL, int, 0, "Number of neutral fitness functions"),
//   VALUE(N_GOOD, int, 7, "Number of good fitness functions"),
//   VALUE(N_BAD, int, 3, "Number of bad fitness functions"),
//   VALUE(DISTANCE_CUTOFF, double, .1, "How close to origin does fitness gradient start"),
  VALUE(RESOURCE_INFLOW, double, 100, "How much resource enters the world each update"),
  VALUE(MUT_RATE, double, .05, "Standard deviation of normal distribution mutations are seelcted from"),
  VALUE(PROBLEM_DIMENSIONS, int, 10, "How many axes does the box have?"),
//   VALUE(RECOMBINATION, int, 0, "Does recombination happen?"),
  VALUE(TOURNAMENT_SIZE, int, 5, "Tournament size"),
  VALUE(COST, double, 0, "Cost of doing task unsuccessfully"),
  VALUE(FRAC, double, .0025, "Percent of resource individual can use"),
  VALUE(GENOME_SIZE, int, 200, "Length of genome"),
  VALUE(MAX_RES_USE, double, 5, "Maximum quantity of resource that individual can use")
)

struct genome_info {
    std::set<size_t> always_used;
    std::set<size_t> sometimes_used;
    bool selected = false;
    emp::vector<double> error_vec;

};

using ORG_TYPE = emp::AvidaGP;

class EcologyWorld  : public emp::World<ORG_TYPE>{
public:

    BoxConfig config;
    uint32_t POP_SIZE;
    uint32_t START_POP_SIZE;
    uint32_t UPDATES;
    double RESOURCE_INFLOW;
    double MUT_RATE;
    std::string SELECTION;
    bool RECOMBINATION;
    int TOURNAMENT_SIZE;
    int GENOME_SIZE;
    double COST;
    double FRAC;
    double MAX_RES_USE;
    TestcaseSet<int, int> testcases;
    std::set<size_t> full_set;
    int failed = 0;

    std::map<ORG_TYPE::genome_t, genome_info> per_genotype_data;


    emp::vector<emp::Resource> resources;
    emp::vector<fun_calc_fitness_t> fit_set;

    emp::Signal<void(emp::vector<size_t>, int repro_id)> on_lexicase_select; ///< Trigger on Lexicase Selection event

    void TriggerOnLexicaseSelect(emp::vector<size_t> vec, int repro_id) {
        // std::cout << "triggering" << std::endl;
        on_lexicase_select.Trigger(vec, repro_id);
    }

    emp::DataNode<size_t, emp::data::Range> niche_width;
    emp::DataNode<size_t, emp::data::Range> sometimes_used;
    emp::DataNode<size_t, emp::data::Range> always_used;

    EcologyWorld(std::string filename) : testcases(filename) {

        for (size_t i = 0; i < testcases.GetTestcases().size(); i++) {
            full_set.insert(i);
        }

        on_lexicase_select.AddAction([this](emp::vector<size_t> vec, int repro_id){
            // std::cout << to_string(vec) << std::endl;
            // std::cout << vec.size() << std::endl;
            niche_width.Add(vec.size());
            // std::cout << niche_width.GetMean() << std::endl;
            per_genotype_data[pop[repro_id]->GetGenome()].sometimes_used.insert(vec.begin(), vec.end());
            per_genotype_data[pop[repro_id]->GetGenome()].always_used = emp::intersection(per_genotype_data[pop[repro_id]->GetGenome()].always_used, vec);
            per_genotype_data[pop[repro_id]->GetGenome()].selected = true;
        });

        OnOffspringReady([this](ORG_TYPE & org) {
            uint32_t num_muts = random_ptr->GetUInt(MUT_RATE*org.GetGenome().sequence.size());  // 0 to 3 mutations.
            for (uint32_t m = 0; m < num_muts; m++) {
                const uint32_t pos = random_ptr->GetUInt(GENOME_SIZE);
                org.RandomizeInst(pos, *random_ptr);
            }
            return num_muts;
        });

        OnOffspringReady([this](ORG_TYPE & org){
            if (!emp::Has(per_genotype_data, org.GetGenome())) {
                per_genotype_data[org.GetGenome()] = genome_info();
                if (SELECTION == "LEXICASE") {
                    per_genotype_data[org.GetGenome()].always_used = full_set;
                }
            }
        });

        OnInjectReady([this](ORG_TYPE & org){
            if (!emp::Has(per_genotype_data, org.GetGenome())) {
                per_genotype_data[org.GetGenome()] = genome_info();
                if (SELECTION == "LEXICASE") {
                    per_genotype_data[org.GetGenome()].always_used = full_set;
                }
            }
        });        
    }

    void SetupFitnessFunctions() {

        fit_set.resize(0);

        int count = 0;
        for (auto testcase : testcases.GetTestcases()) {
            fit_set.push_back([testcase, count, this](ORG_TYPE & org){

                if (emp::Has(per_genotype_data, org.GetGenome()) && (int)per_genotype_data[org.GetGenome()].error_vec.size() > count) {
                    return per_genotype_data[org.GetGenome()].error_vec[count];
                }

                org.ResetHardware();
                for (size_t i = 0; i < testcase.first.size(); i++) {
                    org.SetInput(i, testcase.first[i]);
                }

                org.Process(200);
                int divisor = testcase.second;
                if (divisor == 0) {
                    divisor = 1;
                }   
                double result = 1 - (std::abs(org.GetOutput(0) - testcase.second)/divisor);
                per_genotype_data[org.GetGenome()].error_vec.push_back(result);
                return result;
            });

            count++;
        }


        std::function<double(const ORG_TYPE&)> goal_function = [this](const ORG_TYPE & org){

            double total = 0;
            for (double val : per_genotype_data[org.GetGenome()].error_vec) {
                total += val;
            }

            return total;
        };

        SetFitFun(goal_function);

        // if (SELECTION == "LEXICASE") {
        //     fit_set.push_back(goal_function);
        // } else if (SELECTION == "TOURNAMENT") {
        //     SetCache(true);
        // }

    }


    void InitConfigs() {
        POP_SIZE = config.POP_SIZE();
        START_POP_SIZE = config.START_POP_SIZE();        
        UPDATES = config.UPDATES();
        RESOURCE_INFLOW = config.RESOURCE_INFLOW();
        MUT_RATE = config.MUT_RATE();
        SELECTION = config.SELECTION();
        TOURNAMENT_SIZE = config.TOURNAMENT_SIZE();
        GENOME_SIZE = config.GENOME_SIZE();
        COST = config.COST();
        FRAC = config.FRAC();
        MAX_RES_USE = config.MAX_RES_USE();
    }

    void InitPop() {
        emp::Random & random = GetRandom();
        for (size_t i = 0; i < START_POP_SIZE; i++) {
            ORG_TYPE cpu;
            cpu.PushRandom(random, GENOME_SIZE);
            Inject(cpu.GetGenome());
        }
    }

    void Setup() {

        SetWellMixed(true);

        InitConfigs();

        resources.resize(0);
        for (size_t i=0; i<testcases.GetTestcases().size(); i++) {
            resources.push_back(emp::Resource(RESOURCE_INFLOW, RESOURCE_INFLOW, .01));
        }

        InitPop();

        // Setup the mutation function.
        SetMutFun( [this](ORG_TYPE & org, emp::Random & random) {
            uint32_t num_muts = random.GetUInt(4);  // 0 to 3 mutations.
            for (uint32_t m = 0; m < num_muts; m++) {
                const uint32_t pos = random.GetUInt(GENOME_SIZE);
                org.RandomizeInst(pos, random);
            }
            return num_muts;
        } );

        SetupFitnessFunctions();

        SetupFitnessFile().SetTimingRepeat(10);
        SetupSystematicsFile().SetTimingRepeat(10);
        
        SetupPopulationFile().SetTimingRepeat(10);
        emp::World_file & sel_file = SetupFile("selection_info.csv");
        sel_file.SetTimingRepeat(10);
        sel_file.AddVar(update, "update", "Update");
        sel_file.AddMean(niche_width, "niche_width_mean");
        sel_file.AddMin(niche_width, "niche_width_min");
        sel_file.AddMax(niche_width, "niche_width_max");
        sel_file.AddMean(always_used, "always_used_mean");
        sel_file.AddMin(always_used, "always_used_min");
        sel_file.AddMax(always_used, "always_used_max");
        sel_file.AddMean(sometimes_used, "sometimes_used_mean");
        sel_file.AddMin(sometimes_used, "sometimes_used_min");
        sel_file.AddMax(sometimes_used, "sometimes_used_max");
        sel_file.AddVar(failed, "failed_to_reproduce");
        sel_file.PrintHeaderKeys();

        emp::World_file & extra_systematics_file = SetupFile("extra_systematics.csv");
        extra_systematics_file.SetTimingRepeat(10);
        extra_systematics_file.AddVar(update, "update", "Update");
        extra_systematics_file.AddFun<double>([this](){return systematics.GetPhylogeneticDiversity();}, "phylogenetic_diversity");
        extra_systematics_file.AddFun<double>([this](){return systematics.GetMeanPairwiseDistance();}, "mean_pairwise_distance"); 

        extra_systematics_file.PrintHeaderKeys();

    }

    void RunStep() {
        niche_width.Reset();
        always_used.Reset();
        sometimes_used.Reset();

        EliteSelect(*this);

        if (SELECTION == "TOURNAMENT") {
            TournamentSelect(*this, TOURNAMENT_SIZE, POP_SIZE-1);
        } else if (SELECTION == "LEXICASE") {
            emp::LexicaseSelect(*this, fit_set, POP_SIZE-1);
        } else if (SELECTION == "RESOURCE") {
            emp::ResourceSelect(*this, fit_set, resources, TOURNAMENT_SIZE, POP_SIZE-1, FRAC, MAX_RES_USE, RESOURCE_INFLOW, COST, false);
            for (emp::Ptr<ORG_TYPE> org : pop) {
                if (per_genotype_data[org->GetGenome()].always_used.size() == 0) {
                    std::set<size_t> niches;
                    for (size_t i = 0; i < per_genotype_data[org->GetGenome()].error_vec.size(); i++) {
                        if (per_genotype_data[org->GetGenome()].error_vec[i] > 0) {
                            niches.insert(i);
                        }
                    }
                    per_genotype_data[org->GetGenome()].always_used = niches;
                }
                niche_width.Add(per_genotype_data[org->GetGenome()].always_used.size());
            }
        } else if (SELECTION == "ROULETTE") {
            RouletteSelect(*this, POP_SIZE-1);
        } else {
            std::cout << "ERROR: INVALID SELECTION SCHEME: " << SELECTION << std::endl;
            exit(1);
        }
        
        failed = 0;
        for (emp::Ptr<ORG_TYPE> org : pop) {
            if (!per_genotype_data[org->GetGenome()].selected) {
                failed++;
                continue;
            }
            always_used.Add(per_genotype_data[org->GetGenome()].always_used.size());
            sometimes_used.Add(per_genotype_data[org->GetGenome()].sometimes_used.size());
        }
        // std::cout << niche_width.GetMax() << std::endl;
        
        std::cout << "UD: " << update << " Best: " << GetFitnessDataNode().GetMax() << std::endl;

        Update();
        // if (isinf(GetFitnessDataNode().GetMax())){
        //     return;
        // }
        
        // for (auto res : resources) {
        //     std::cout << res.GetAmount() << " ";
        // }
        // std::cout << std::endl;
    }

    void Run() {
        for (size_t ud = 0; ud < UPDATES; ud++) {
            RunStep();
            // if (isinf(GetFitnessDataNode().GetMax())){
            //     return;
            // }
        }
        PrintDetail("detail.spop");
    }

    void PrintDetail(std::string filename) {
        std::ofstream detail_file;
        detail_file.open(filename);

        if (!detail_file.is_open()) {
            std::cout << "ERROR: failed to open detail file " << filename << std::endl;
        }

        for (const auto & x : systematics.GetActive()) {
            int parent_id;
            if (x->GetParent()) {
                parent_id = x->GetParent();
            }

          detail_file << x->GetID() << "," << x->GetNumOrgs() << "," << x->GetNumOff() << ","
             << parent_id << "," << emp::to_string(per_genotype_data[x->GetInfo()].error_vec) <<
             "," << emp::to_string(per_genotype_data[x->GetInfo()].always_used) << "," << 
             emp::to_string(per_genotype_data[x->GetInfo()].sometimes_used) << std::endl;
        }


        for (const auto & x : systematics.GetAncestors()) {
            int parent_id;
            if (x->GetParent()) {
                parent_id = x->GetParent();
            }
            detail_file << x->GetID() << "," << x->GetNumOrgs() << "," << x->GetNumOff() << "|"
             << parent_id << "," << emp::to_string(per_genotype_data[x->GetInfo()].error_vec) <<
             "," << emp::to_string(per_genotype_data[x->GetInfo()].always_used) << "," << 
             emp::to_string(per_genotype_data[x->GetInfo()].sometimes_used) << std::endl;
        }



    }
};
