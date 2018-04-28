#include "config/ArgManager.h"
#include "base/vector.h"
#include "tools/Random.h"
#include "tools/math.h"
#include "tools/set_utils.h"
#include "tools/stats.h"
#include "Evolve/World.h"
#include "Evolve/Resource.h"
#include "hardware/AvidaGP.h"
#include "TestcaseSet.h"
#include "data/DataNode.h"
#include "control/Signal.h"
#include "tools/map_utils.h"
#include "tools/memo_function.h"
#include <map>
#include <unordered_map>

// #include "cec2013.h"

namespace std {
    template <>
    struct hash<std::vector<double>> {
        std::size_t operator() (const std::vector<double> &vc) const
        {
            std::size_t seed = vc.size();
            for (auto& i : vc) {
                seed ^= std::hash<double>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    #ifndef NDEBUG

    template <>
    struct hash<emp::vector<double>> {
        std::size_t operator() (const emp::vector<double> &vc) const
        {
            std::size_t seed = vc.size();
            for (auto& i : vc) {
                seed ^= std::hash<double>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    #endif

    template <>
    struct hash<emp::AvidaCPU_Base<emp::AvidaGP>::Genome> {
        std::size_t operator() (const emp::AvidaCPU_Base<emp::AvidaGP>::Genome &vc) const
        {
            return vc.Hash();
        }
    };
}

EMP_BUILD_CONFIG( BoxConfig,
  GROUP(DEFAULT, "Default settings for box experiment"),
  VALUE(SEED, int, 55711224, "Random number seed (0 for based on time)"),
  VALUE(POP_SIZE, uint32_t, 1000, "Number of organisms in the popoulation."),
  VALUE(START_POP_SIZE, uint32_t, 1, "Number of organisms to initialize population with."),
  VALUE(UPDATES, uint32_t, 1000, "How many generations should we process?"),
  VALUE(SELECTION, std::string, "LEXICASE", "What selection scheme should we use?"),
  VALUE(N_NEUTRAL, int, 0, "Number of neutral fitness functions"),
  VALUE(N_GOOD, int, 7, "Number of good fitness functions"),
  VALUE(N_BAD, int, 3, "Number of bad fitness functions"),
  VALUE(DISTANCE_CUTOFF, double, .1, "How close to origin does fitness gradient start"),
  VALUE(RESOURCE_INFLOW, double, 100, "How much resource enters the world each update"),
  VALUE(MUT_RATE, double, .02, "Standard deviation of normal distribution mutations are seelcted from"),
  VALUE(PROBLEM_DIMENSIONS, int, 10, "How many axes does the box have?"),
  VALUE(PROBLEM, std::string, "testcases/count-odds.csv", "Which set of testcases should we use? (or enter 'box' for the box problem"),
//   VALUE(RECOMBINATION, int, 0, "Does recombination happen?"),
  VALUE(TOURNAMENT_SIZE, int, 5, "Tournament size"),
  VALUE(N_TEST_CASES, uint32_t, 200, "How many test cases to use"),  
  VALUE(COST, double, 0, "Cost of doing task unsuccessfully"),
  VALUE(FRAC, double, .0025, "Percent of resource individual can use"),
  VALUE(GENOME_SIZE, int, 200, "Length of genome"),
  VALUE(SHARE_THRESHOLD, double, 10, "How similar do two solutions need to be to be in the same niche?"),  
  VALUE(MAX_RES_USE, double, 5, "Maximum quantity of resource that individual can use")
)

struct genome_info {
    std::set<size_t> always_used;
    std::set<size_t> sometimes_used;
    bool selected = false;
    emp::vector<double> error_vec;
    double fit = 0;

};

template <typename ORG_TYPE>
class EcologyWorld : public emp::World<ORG_TYPE>{
public:

    using typename emp::World<ORG_TYPE>::genotype_t;
    using typename emp::World<ORG_TYPE>::genome_t;
    using typename emp::World<ORG_TYPE>::fun_calc_fitness_t;
    using emp::World<ORG_TYPE>::random_ptr; 
    using emp::World<ORG_TYPE>::GetRandom;
    using emp::World<ORG_TYPE>::SetMutFun;
    using emp::World<ORG_TYPE>::OnInjectReady;
    using emp::World<ORG_TYPE>::OnOffspringReady;
    using emp::World<ORG_TYPE>::update;
    using emp::World<ORG_TYPE>::Update;
    using emp::World<ORG_TYPE>::pop;
    using emp::World<ORG_TYPE>::systematics;
    using emp::World<ORG_TYPE>::SetupFile;
    using emp::World<ORG_TYPE>::SetupFitnessFile;
    using emp::World<ORG_TYPE>::SetupSystematicsFile;
    using emp::World<ORG_TYPE>::SetupPopulationFile;
    using emp::World<ORG_TYPE>::GetFitnessDataNode;    
    using emp::World<ORG_TYPE>::SetWellMixed;
    using emp::World<ORG_TYPE>::Inject;    
    using emp::World<ORG_TYPE>::GetGenomeAt;    
    using emp::World<ORG_TYPE>::GetGenome;    
    using emp::World<ORG_TYPE>::SetCache;    
    using emp::World<ORG_TYPE>::ClearCache;    

    uint32_t POP_SIZE;
    uint32_t START_POP_SIZE;
    int N_NEUTRAL;
    int N_GOOD;
    int N_BAD;
    int PROBLEM_DIMENSIONS;
    double DISTANCE_CUTOFF;
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
    TestcaseSet<int, double> testcases;
    std::set<size_t> full_set;
    int failed = 0;
    std::string PROBLEM;
    uint32_t N_TEST_CASES;
    double SHARE_THRESHOLD;

    bool do_mutate = true;

    std::unordered_map<genome_t, genome_info> per_genotype_data;

    bool destructed = false;

    emp::vector<emp::Resource> resources;
    emp::vector<fun_calc_fitness_t> fit_set;

    emp::Signal<void(emp::vector<size_t>, int repro_id)> on_lexicase_select; ///< Trigger on Lexicase Selection event
    // emp::Signal<void(int repro_id)> on_elite_select; ///< Trigger on Lexicase Selection event

    void TriggerOnLexicaseSelect(emp::vector<size_t> vec, int repro_id) {
        // std::cout << "triggering" << std::endl;
        on_lexicase_select.Trigger(vec, repro_id);
    }

    // void TriggerOnEliteSelect(int repro_id) {
    //     // std::cout << "triggering" << std::endl;
    //     on_elite_select.Trigger(repro_id);
    // }


    emp::DataNode<size_t, emp::data::Range> niche_width;
    emp::DataNode<size_t, emp::data::Range> sometimes_used;
    emp::DataNode<size_t, emp::data::Range> always_used;

    emp::DataNode<double, emp::data::Range> evolutionary_distinctiveness;


    EcologyWorld() {}
    EcologyWorld(emp::Random & rnd) : emp::World<ORG_TYPE>(rnd) {}
    ~EcologyWorld() {destructed = true;};



    void InitConfigs(BoxConfig & config) {
        POP_SIZE = config.POP_SIZE();
        START_POP_SIZE = config.START_POP_SIZE();        
        UPDATES = config.UPDATES();
        N_NEUTRAL = config.N_NEUTRAL();
        N_GOOD = config.N_GOOD();
        N_BAD = config.N_BAD();
        PROBLEM_DIMENSIONS = config.PROBLEM_DIMENSIONS();
        DISTANCE_CUTOFF = config.DISTANCE_CUTOFF();        
        RESOURCE_INFLOW = config.RESOURCE_INFLOW();
        MUT_RATE = config.MUT_RATE();
        SELECTION = config.SELECTION();
        TOURNAMENT_SIZE = config.TOURNAMENT_SIZE();
        GENOME_SIZE = config.GENOME_SIZE();
        COST = config.COST();
        FRAC = config.FRAC();
        MAX_RES_USE = config.MAX_RES_USE();
        PROBLEM = config.PROBLEM();
        N_TEST_CASES = config.N_TEST_CASES();        
        SHARE_THRESHOLD = config.SHARE_THRESHOLD();        
    }
    void InitPop();

    void SetupFitnessFunctions();

    void SetupMutationFunctions();
    void PrintBest();

    double PhenotypicRichness() {
        std::map<emp::vector<double>, int> phen_counts;
        for (auto tax : systematics.GetActive()) {
            emp::vector<double> phen = per_genotype_data[tax->GetInfo()].error_vec;
            if (emp::Has(phen_counts, phen)) {
                phen_counts[phen] += tax->GetNumOrgs();
            } else {
                phen_counts[phen] = tax->GetNumOrgs();
            }
        }

        return phen_counts.size();
    }

    double PhenotypicEntropy() {
        std::map<emp::vector<double>, int> phen_counts;
        for (auto tax : systematics.GetActive()) {
            emp::vector<double> phen = per_genotype_data[tax->GetInfo()].error_vec;
            if (emp::Has(phen_counts, phen)) {
                phen_counts[phen] += tax->GetNumOrgs();
            } else {
                phen_counts[phen] = tax->GetNumOrgs();
            }
        }

        emp::vector<int> counts;
        for (auto phen : phen_counts) {
            counts.push_back(phen.second);
        }

        return emp::Entropy(counts);
    }

    void Setup(BoxConfig & config) {

        SetWellMixed(true);
        SetCache(true);

        InitConfigs(config);

        SetupFitnessFunctions();

        for (size_t i = 0; i < N_TEST_CASES; i++) {
            full_set.insert(i);
        }

        on_lexicase_select.AddAction([this](emp::vector<size_t> vec, int repro_id){
            // std::cout << to_string(vec) << std::endl;
            // std::cout << vec.size() << std::endl;
            niche_width.Add(vec.size());
            // std::cout << niche_width.GetMean() << std::endl;
            per_genotype_data[GetGenomeAt(repro_id)].sometimes_used.insert(vec.begin(), vec.end());
            per_genotype_data[GetGenomeAt(repro_id)].always_used = emp::intersection(per_genotype_data[GetGenomeAt(repro_id)].always_used, vec);
            per_genotype_data[GetGenomeAt(repro_id)].selected = true;
        });

        SetupMutationFunctions();

        resources.resize(0);
        for (size_t i=0; i<N_TEST_CASES; i++) {
            resources.push_back(emp::Resource(RESOURCE_INFLOW, RESOURCE_INFLOW, .01));
        }

        systematics.OnNew([this](emp::Ptr<genotype_t> org){
            if (!emp::Has(per_genotype_data, org->GetInfo())){

                genome_info gen;
                ORG_TYPE cpu = org->GetInfo();
                for (auto fit_fun : fit_set) {
                    gen.error_vec.push_back(fit_fun(cpu));
                }
                gen.fit = emp::Sum(gen.error_vec);

                if (SELECTION == "LEXICASE") {
                    gen.always_used = full_set;
                }
                per_genotype_data[org->GetInfo()] = gen;
            }

        });

        systematics.OnPrune([this](emp::Ptr<genotype_t> org){
            
            if (!destructed && emp::Has(per_genotype_data, org->GetInfo())) {
                per_genotype_data.erase(org->GetInfo());
            }

        });


        InitPop();

        // // Setup the mutation function.
        // SetMutFun( [this](ORG_TYPE & org, emp::Random & random) {
        //     if (!do_mutate) {
        //         return (uint32_t)0;
        //     }
        //     uint32_t num_muts = random.GetUInt(4);  // 0 to 3 mutations.
        //     for (uint32_t m = 0; m < num_muts; m++) {
        //         const uint32_t pos = random.GetUInt(GENOME_SIZE);
        //         org.RandomizeInst(pos, random);
        //     }
        //     return num_muts;
        // } );

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
        extra_systematics_file.AddMean(evolutionary_distinctiveness, "mean_evolutionary_distinctiveness"); 
        extra_systematics_file.AddMin(evolutionary_distinctiveness, "min_evolutionary_distinctiveness");     
        extra_systematics_file.AddMax(evolutionary_distinctiveness, "max_evolutionary_distinctiveness"); 
        extra_systematics_file.AddFun<double>([this](){return PhenotypicRichness();}, "phenotypic_richness"); 
        extra_systematics_file.AddFun<double>([this](){return PhenotypicEntropy();}, "phenotypic_entropy");         
        
        extra_systematics_file.PrintHeaderKeys();

    }

    void RunStep() {
        niche_width.Reset();
        always_used.Reset();
        sometimes_used.Reset();
        evolutionary_distinctiveness.Reset();

        do_mutate = false;
        EliteSelect(*this);
        do_mutate = true;

        if (SELECTION == "TOURNAMENT" || SELECTION == "SHARING") {
            TournamentSelect(*this, TOURNAMENT_SIZE, POP_SIZE-1);
        } else if (SELECTION == "LEXICASE") {
	  std::cout << "lex"<<std::endl;
            emp::LexicaseSelect(*this, fit_set, POP_SIZE-1);
        } else if (SELECTION == "RESOURCE") {
            emp::ResourceSelect(*this, fit_set, resources, TOURNAMENT_SIZE, POP_SIZE-1, FRAC, MAX_RES_USE, RESOURCE_INFLOW, COST, true);
            for (emp::Ptr<ORG_TYPE> org : pop) {
                if (per_genotype_data[GetGenome(*org)].always_used.size() == 0) {
                    std::set<size_t> niches;
                    for (size_t i = 0; i < per_genotype_data[GetGenome(*org)].error_vec.size(); i++) {
                        if (per_genotype_data[GetGenome(*org)].error_vec[i] > 0) {
                            niches.insert(i);
                        }
                    }
                    per_genotype_data[GetGenome(*org)].always_used = niches;
                }
                niche_width.Add(per_genotype_data[GetGenome(*org)].always_used.size());
            }
        } else if (SELECTION == "ROULETTE") {
            RouletteSelect(*this, POP_SIZE-1);
        } else {
            std::cout << "ERROR: INVALID SELECTION SCHEME: " << SELECTION << std::endl;
            exit(1);
        }
        
        failed = 0;
        for (emp::Ptr<ORG_TYPE> org : pop) {
            if (!per_genotype_data[GetGenome(*org)].selected) {
                failed++;
                continue;
            }
            always_used.Add(per_genotype_data[GetGenome(*org)].always_used.size());
            sometimes_used.Add(per_genotype_data[GetGenome(*org)].sometimes_used.size());
        }

        for (auto tax : systematics.GetActive() ) {
            evolutionary_distinctiveness.Add(systematics.GetEvolutionaryDistinctiveness(tax, update));
        }
        // std::cout << niche_width.GetMax() << std::endl;
        
        std::cout << "UD: " << update << " | Best: " << GetFitnessDataNode().GetMax() << std::endl;

        if (update % 100 == 0) {
            PrintDetail("detail-" + emp::to_string(update) +".spop");
        }

        Update();
        ClearCache();
        PrintBest();
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
    }

    void PrintDetail(std::string filename) {
        std::ofstream detail_file;
        detail_file.open(filename);

        if (!detail_file.is_open()) {
            std::cout << "ERROR: failed to open detail file " << filename << std::endl;
        }

        for (const auto & x : systematics.GetActive()) {
            int parent_id = 0;
            if (x->GetParent()) {
                parent_id = x->GetParent()->GetID();
            }

          detail_file << x->GetID() << "," << x->GetNumOrgs() << "," << x->GetNumOff() << ","
             << parent_id << "," << emp::to_string(per_genotype_data[x->GetInfo()].error_vec) <<
             "," << emp::to_string(per_genotype_data[x->GetInfo()].always_used) << "," << 
             emp::to_string(per_genotype_data[x->GetInfo()].sometimes_used) << std::endl;
        }


        for (const auto & x : systematics.GetAncestors()) {
            int parent_id = 0;
            if (x->GetParent()) {
                parent_id = x->GetParent();
            }
            detail_file << x->GetID() << "," << x->GetNumOrgs() << "," << x->GetNumOff() << ","
             << parent_id << "," << emp::to_string(per_genotype_data[x->GetInfo()].error_vec) <<
             "," << emp::to_string(per_genotype_data[x->GetInfo()].always_used) << "," << 
             emp::to_string(per_genotype_data[x->GetInfo()].sometimes_used) << std::endl;
        }



    }
};

template <> 
void EcologyWorld<emp::AvidaGP>::PrintBest() {
        std::ofstream detail_file;
        detail_file.open("best.org");

        for (auto org : pop) {
            if (emp::Sum(per_genotype_data[org->GetGenome()].error_vec) == per_genotype_data[org->GetGenome()].error_vec.size()) {
                pop[0]->PrintGenome(detail_file);
                detail_file.close();
                return;
   
            }
        }

}

template <> 
void EcologyWorld<emp::vector<double>>::PrintBest() {
        std::ofstream detail_file;
        detail_file.open("best.org");

        detail_file << emp::to_string(*pop[0]);
        detail_file.close();
}

template <typename ORG_TYPE> 
void EcologyWorld<ORG_TYPE>::PrintBest() {
    /* This is an org_type that we don't have any fitness functions for*/
    std::cout << "Unsporrted org type." << std::endl;
    exit(1);
}

template <typename ORG_TYPE> 
void EcologyWorld<ORG_TYPE>::SetupMutationFunctions() {
    /* This is an org_type that we don't have any fitness functions for*/
    std::cout << "Unsporrted org type." << std::endl;
    exit(1);
}


template <> 
void EcologyWorld<emp::AvidaGP>::SetupMutationFunctions() {

    OnOffspringReady([this](emp::AvidaGP & org) {
        if (do_mutate) {
            
            uint32_t num_muts = random_ptr->GetUInt(MUT_RATE);  // 0 to 3 mutations.
            for (uint32_t m = 0; m < num_muts; m++) {
                const uint32_t pos = random_ptr->GetUInt(GENOME_SIZE);
                org.RandomizeInst(pos, *random_ptr);
            }
        }

    });
}


template <> 
void EcologyWorld<emp::vector<double> >::SetupMutationFunctions() {

    OnOffspringReady([this](emp::vector<double> & org) {
        if (do_mutate) {
            for (int pos = 0; pos < GENOME_SIZE; pos++) {
                org[pos] += random_ptr->GetRandNormal(0, MUT_RATE);
                if (org[pos] < 0) {
                    org[pos] = 0;
                } else if (org[pos] > 1) {
                    org[pos] = 1;
                }
            }
        }
    });
}


template <typename ORG_TYPE> 
void EcologyWorld<ORG_TYPE>::InitPop() {
    /* This is an org_type that we don't have any fitness functions for*/
    std::cout << "Unsporrted org type." << std::endl;
    exit(1);
}


template <>
void EcologyWorld<emp::vector<double> >::InitPop() {
    emp::Random & random = GetRandom();
    for (size_t i = 0; i < START_POP_SIZE; i++) {
        emp::vector<double> vec;
        for (int j = 0; j < GENOME_SIZE; j++ ) {
            vec.push_back(random.GetDouble(0, 1));
        }
        Inject(vec);
    }
}

template <>
void EcologyWorld<emp::AvidaGP>::InitPop() {
    emp::Random & random = GetRandom();
    for (size_t i = 0; i < START_POP_SIZE; i++) {
        emp::AvidaGP cpu;
        cpu.PushRandom(random, GENOME_SIZE);
        Inject(cpu.GetGenome());
    }
}


template <typename ORG_TYPE> 
void EcologyWorld<ORG_TYPE>::SetupFitnessFunctions() {
    /* This is an org_type that we don't have any fitness functions for*/
    std::cout << "Unsporrted org type." << std::endl;
    exit(1);
}

template <>
void EcologyWorld<emp::vector<double> >::SetupFitnessFunctions() {

    fit_set.resize(0);

    if (PROBLEM == "box") {
        // Good hints
        for (int i = 0; i < N_GOOD; i++) {
            fit_set.push_back([i](const emp::vector<double> & org){
                double score = 1 - org[i];
                if (score < .8) {
                    return 0.0;
                }
                return score;
            });
        }

        // Bad hints
        for (int i = N_GOOD; i < N_GOOD + N_BAD; i++) {
            fit_set.push_back([i](const emp::vector<double> & org){
                double score = org[i];
                if (score < .8) {
                    return 0.0;
                }
                return score;
            });
        }

        // Neutral hints (these axes aren't evaluated)
        for (int i = PROBLEM_DIMENSIONS; i < GENOME_SIZE; i++) {
            fit_set.push_back([i](const emp::vector<double> & org){return org[i];});
        }

        fun_calc_fitness_t goal_function = [this](const emp::vector<double> & org){
            double dist = 0;
            for (int i = 0; i < PROBLEM_DIMENSIONS; i++) {
                dist += emp::Pow(org[i], 2.0);
            }

            dist = sqrt(dist);
            if (dist > DISTANCE_CUTOFF) {
                return 0.01;
            }

            return 1.0/dist;
        };

        fun_calc_dist_t dist_fun = [this](const emp::vector<double> & org1, const emp::vector<double> & org2) {
            double dist = 0;

            emp_assert(org1.size() == org2.size(), org1.size(), org2.size());

            for (size_t i = 0; i < org1.size(); i++) {
                dist += emp::Pow(org1[i]-org2[i], 2.0);
            }

            return sqrt(dist);
        };
        


        if (SELECTION == "SHARING") {
            SetSharedFitFun(goal_function, dist_fun, .1, 1);
        } else {
            SetFitFun(goal_function);
        }

    } else {
        std::cout << "Invalid problem specified" << std::endl;
    }
    // if (SELECTION == "LEXICASE") {
    //     fit_set.push_back(goal_function);
    // } else if (SELECTION == "TOURNAMENT") {
    //     SetCache(true);
    // }
}

template <>
void EcologyWorld<emp::AvidaGP>::SetupFitnessFunctions() {

    testcases.LoadTestcases(PROBLEM);
    fit_set.resize(0);

    int count = 0;
    for (auto testcase : testcases.GetTestcases()) {
        fit_set.push_back([testcase, count, this](emp::AvidaGP & org){

            if (emp::Has(per_genotype_data, org.GetGenome()) && ((uint32_t)(per_genotype_data[org.GetGenome()].error_vec.size()) == N_TEST_CASES)) {
                return per_genotype_data[org.GetGenome()].error_vec[count];
            }

            org.ResetHardware();
            for (size_t i = 0; i < testcase.first.size(); i++) {
                org.SetInput(i, testcase.first[i]);
            }
            org.SetOutput(0,-99999); // Otherwise not outputting anything is a decent strategy

            org.Process(200);
            int divisor = testcase.second;
            if (divisor == 0) {
                divisor = 1;
            }   
            double result = 1 - (std::abs(org.GetOutput(0) - testcase.second)/divisor);
            // emp_assert(std::abs(result) != INFINITY);
            if (result == -INFINITY) {
                result = -999999999;
            }
            // per_genotype_data[org.GetGenome()].error_vec.push_back(result);
            return result;
        });

        count++;
        if (count >= (int)N_TEST_CASES) {
            break;
        }
    }


    std::function<double(const emp::AvidaGP&)> goal_function = [this](const emp::AvidaGP & org){

        // double total = 0;
        // for (double val : per_genotype_data[org.GetGenome()].error_vec) {
        //     total += val;
        // }

        // return total;
        return per_genotype_data[org.GetGenome()].fit;

    };

    fun_calc_dist_t dist_fun = [this](emp::AvidaGP & org1, emp::AvidaGP & org2) {
        double dist = 0;

        // emp_assert(org1.size() == org2.size(), org1.size(), org2.size());

        // emp_assert(emp::Has(per_genotype_data, GetGenome(org1)));
        // emp_assert(emp::Has(per_genotype_data, GetGenome(org2)));

        auto it1 = per_genotype_data.find(GetGenome(org1));
        auto it2 = per_genotype_data.find(GetGenome(org2));

        if (it1 == per_genotype_data.end()){

            per_genotype_data[GetGenome(org1)] = genome_info();
            for (auto fit_fun : fit_set) {
                per_genotype_data[GetGenome(org1)].error_vec.push_back(fit_fun(org1));
            }

            per_genotype_data[GetGenome(org1)].fit = emp::Sum(per_genotype_data[GetGenome(org1)].error_vec);

            if (SELECTION == "LEXICASE") {
                per_genotype_data[GetGenome(org1)].always_used = full_set;
            }
            it1 = per_genotype_data.find(GetGenome(org1));
        }

        if (it2 == per_genotype_data.end()){

            per_genotype_data[GetGenome(org2)] = genome_info();
            for (auto fit_fun : fit_set) {
                per_genotype_data[GetGenome(org2)].error_vec.push_back(fit_fun(org2));
            }
            per_genotype_data[GetGenome(org1)].fit = emp::Sum(per_genotype_data[GetGenome(org1)].error_vec);

            if (SELECTION == "LEXICASE") {
                per_genotype_data[GetGenome(org2)].always_used = full_set;
            }
            it2 = per_genotype_data.find(GetGenome(org2));
        }


        emp::vector<double> & err_vec1 = it1->second.error_vec;
        emp::vector<double> & err_vec2 = it2->second.error_vec;

        for (size_t i = 0; i < N_TEST_CASES; i++) {
            dist += emp::Pow(err_vec1[i] - err_vec2[i], 2.0);
        }

        return sqrt(dist);
    };
    


    if (SELECTION == "SHARING") {
        SetSharedFitFun(goal_function, dist_fun, SHARE_THRESHOLD, 1);
    } else {
        SetFitFun(goal_function);
    }
}
