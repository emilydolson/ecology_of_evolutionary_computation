#ifndef _INTERACTION_NETWORKS_H_
#define _INTERACTION_NETWORKS_H_

#include <functional>
#include <map>
#include <limits>

#include "web/init.h"

#include "base/vector.h"
#include "base/assert.h"
#include "tools/Random.h"
#include "tools/random_utils.h"
#include "tools/vector_utils.h"
#include "tools/map_utils.h"
#include "tools/string_utils.h"
#include "tools/math.h"
#include "tools/distances.h"
#include "tools/attrs.h"
#include "tools/Graph.h"
#include "Evolve/Resource.h"

#include <permutations.hpp>

using org_t = emp::vector<int>;
using fit_map_t = std::map<org_t, double>;


DEFINE_ATTR(SigmaShare);
DEFINE_ATTR(Alpha);
DEFINE_ATTR(Cost);
DEFINE_ATTR(Cf);
DEFINE_ATTR(NicheWidth);
DEFINE_ATTR(MaxScore);
DEFINE_ATTR(ResourceInflow);
DEFINE_ATTR(ResourceOutflow);
DEFINE_ATTR(MaxBonus);

constexpr auto DEFAULT{MakeAttrs(SigmaShare(8.0),
                                 Alpha(1.0),
                                 Cost(1.0),
                                 Cf(.0025),                                 
                                 NicheWidth(3.0),
                                 MaxScore(10.0),
                                 ResourceInflow(2000),
                                 ResourceOutflow(.01),
                                 MaxBonus(5))};

using all_attrs = emp::tools::Attrs<typename SigmaShare::value_t<double>, typename Alpha::value_t<double>, 
                        typename Cost::value_t<double>, typename Cf::value_t<double>,
                        typename NicheWidth::value_t<double>, typename MaxScore::value_t<double>,
                        typename ResourceInflow::value_t<double>, typename ResourceOutflow::value_t<double>,
                        typename MaxBonus::value_t<double> >;


// from https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(T x, T y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
    // unless the result is subnormal
           || std::abs(x-y) < std::numeric_limits<T>::min();
}

emp::vector<org_t> make_pop(emp::Random & r, 
                                        int size = 10, 
                                        int length=5, 
                                        double p=.5) {
    emp::vector< org_t> pop;
    pop.resize(size);

    for (int org = 0; org < size; ++org) {
        for (int loc = 0; loc < length; ++loc) {
            pop[org].push_back(r.GetRandGeometric(p));
        }
    } 
    return pop;
}

emp::vector<org_t> FindHighest(emp::vector<org_t> & pop, int axis) {
    double best = -1;
    emp::vector<org_t> winners;

    for (org_t & org : pop) {
        if (org[axis] > best) {
            best = org[axis];
            winners.resize(0);
            winners.push_back(org);
        } else if (org[axis] == best) {
            winners.push_back(org);
        }
    }
    return winners;
}

long int PartialFactorial(int start, int stop) {
    // std::cout << start << " " << stop << std::endl;
    long int result = start;
    while (start > stop) {
        result *= start - 1;
        start--;        
        // std::cout << "test: " << start << std::endl;
    }
    // std::cout << result << std::endl;
    return result;
}

void TraverseDecisionTree(fit_map_t & fit_map, emp::vector<org_t> & pop, emp::vector<int> axes) {
    if (axes.size() == 1) {
        emp::vector<org_t> winners = FindHighest(pop, axes[0]);
        for (org_t & org : winners) {
            fit_map[org]+=1.0/((double)winners.size()*emp::Factorial(pop[0].size()));
        }
        return;
    }

    for (int ax : axes) {
        // std::cout << "Axis: " << ax << " out of " << emp::to_string(axes) << std::endl;
        emp::vector<org_t> winners = FindHighest(pop, ax);
        if (winners.size() == 1) { // Not a tie
            // std::cout << "1 winner: " << emp::to_string(winners[0]) << " Controls " << (double)emp::Factorial(axes.size() - 1)<< std::endl;
            fit_map[winners[0]] += 1.0/PartialFactorial(pop[0].size(), axes.size());//(double)emp::Factorial(axes.size() - 1);
        } else { // tie
            emp::vector<int> next_axes;
            for (int new_ax : axes) {
                if (new_ax != ax) {
                    next_axes.push_back(new_ax);
                }
            }
            TraverseDecisionTree(fit_map, winners, next_axes);
        }
    }
}


std::function<fit_map_t(emp::vector<org_t>&, all_attrs)> lexicase_fitness = [](emp::vector<org_t> & pop, all_attrs attrs=DEFAULT) {

    // std::cout << "LEXICASE" << std::endl;
    emp_assert(pop.size() > 0);
    fit_map_t fit_map;
    size_t n_funs = pop[0].size();

    for (org_t & org : pop) {
        fit_map[org] = 0.0;
    }

    emp::vector<org_t> de_dup_pop = emp::RemoveDuplicates(pop);
    TraverseDecisionTree(fit_map, de_dup_pop, emp::NRange(0, (int)n_funs));

    for (org_t & org : de_dup_pop) {
        fit_map[org] /= emp::Count(pop, org);
        // std::cout << "Pre div: " << fit_map[org] << std::endl;
        // fit_map[org] /= emp::Factorial(n_funs); // convert to proportion of "islands"
        // std::cout << "Post div: " << fit_map[org] << " Divided by: " << (double)emp::Factorial(n_funs) << std::endl;
    }

    return fit_map;
};

std::function<fit_map_t(emp::vector<org_t>&, all_attrs)> sharing_fitness = [](emp::vector<org_t> & pop, all_attrs attrs=DEFAULT) {
    // std::cout << "SHARING" << std::endl;
    fit_map_t base_fit_map;
    fit_map_t fit_map;

    // std::cout << SigmaShare::Get(attrs) << std::endl;

    for (org_t & org : pop) {
        base_fit_map[org] = 1.0;
        fit_map[org] = 0.0;
    }

    for (org_t & org1 : pop) {
        double niche_count = 0;

        for (org_t & org2 : pop) {
            // Sharing function is euclidean distance
            // we could make this configurable
            double dist = emp::EuclideanDistance(org1, org2);
            if (dist < SigmaShare::Get(attrs)) {
                niche_count += 1 - emp::Pow((dist/SigmaShare::Get(attrs)), Alpha::Get(attrs));
            } 
        }

        // Don't worry, niche_count will always be at least one because each individual
        // increases its own niche count by 1
        base_fit_map[org1] = emp::Sum(org1) / niche_count;
    }

    // Now that we've adjusted for sharing, calculate biological fitness
    for (org_t & org : pop) {
        double less = 0.0;
        double equal = -1; // Org will be equal to itself
        double greater = 0.0;

        for (org_t & org2 : pop) {
            if (base_fit_map[org2] < base_fit_map[org]) {
                less++;
            } else if (almost_equal(base_fit_map[org2], base_fit_map[org], 10)) {
                equal++;
            } else {
                greater++;
            }
        }

        double p_less = less/((int)pop.size() - 1);
        double p_equal = equal/((int)pop.size() - 1);
        // double p_greater = greater/(int)pop.size();
        
        fit_map[org] = (2*p_less + p_equal) / (int)pop.size();
    }

    return fit_map;    
};


std::function<fit_map_t(emp::vector<org_t>&, all_attrs)> eco_ea_fitness = [](emp::vector<org_t> & pop, all_attrs attrs=DEFAULT) {
    // std::cout << "ECO-EA" << std::endl;
    emp_assert(pop.size() > 0);

    fit_map_t base_fit_map;
    fit_map_t fit_map;

    size_t n_funs = pop[0].size();

    for (org_t & org : pop) {
        base_fit_map[org] = 1.0;
        fit_map[org] = 0.0;
    }

    for (int axis : emp::NRange(0, (int)n_funs)) {
        double res = ResourceInflow::Get(attrs);
        double count = 0;

        for (org_t & org : pop) {
            if (org[axis] >= NicheWidth::Get(attrs)) {
                count++;
            }
        }

        if (count > 0) {
            res /= count; // Ignores resource accumulation, but that's okay for interactions
        }

        for (org_t & org : pop) {
            if (org[axis] >= NicheWidth::Get(attrs)) {
                base_fit_map[org] *= emp::Pow2(std::min(Cf::Get(attrs)*res*emp::Pow(org[axis]/MaxScore::Get(attrs),2.0) - Cost::Get(attrs), MaxBonus::Get(attrs)));
            }
        }
    }

    for (org_t & org : emp::RemoveDuplicates(pop)) {
        int wins = 0;
        int ties = -1; // org will tie with itself
        double fit = base_fit_map[org];

        for (org_t & org2 : pop) {
            // std::cout << fit << ' ' << base_fit_map[org2] << std::endl;
            if (fit > base_fit_map[org2]) {
                // std::cout << "win" << std::endl;
                wins++;
            } else if (almost_equal(fit, base_fit_map[org2], 10)) {
                // std::cout << "tie " << ties << std::endl;
                ties++;
                // std::cout << ties << std::endl;
            }
        }

        // Assumes tournament size of 2
        fit_map[org] = (2/(double)pop.size()) * (wins/((double)pop.size() - 1) + .5*ties/((double)pop.size() - 1));
        // std::cout << "IN ECOEA " << fit_map[org] << " wins: " << wins << " ties: " << ties << std::endl;
    }

    return fit_map;
};

emp::WeightedGraph calc_competition(emp::vector<org_t> pop, 
                        std::function<fit_map_t(emp::vector<org_t>&, all_attrs)> fit_fun,
                        all_attrs attrs=DEFAULT) {

    emp::WeightedGraph effects(pop.size());


    fit_map_t fitnesses = fit_fun(pop, attrs);

    for (size_t i = 0; i < pop.size(); i++) {
        effects.SetLabel(i, emp::to_string(pop[i]));
        // std::cout << effects.GetLabel(i) << std::endl;

        emp::vector<org_t> curr = pop;
        for (int & ax : curr[i]) {
            ax = 0; // Replace org with null org so pop size stays same
        }

        fit_map_t new_fits = fit_fun(curr, attrs);
        for (size_t j = 0; j < pop.size(); j++ ) {
            if (i == j) {continue;}
            double effect = fitnesses[curr[j]] - new_fits[curr[j]];
            if (effect != 0) {
                effects.AddEdge(i, j, effect);
            }
            // std::cout << effect << std::endl;
        }
    }

    return effects;
}

class Controller {
public:

    emp::vector<org_t> pop;
    int pop_size = 10;
    int n_traits = 5;
    all_attrs settings = DEFAULT;
    emp::Random r;
    emp::vector<emp::Resource> resources;

    emp::WeightedGraph eco_network;
    emp::WeightedGraph lex_network;
    emp::WeightedGraph share_network;
    
    Controller(int pop=10, int traits=5) : pop_size(pop), n_traits(traits) {
        EM_ASM((window["js"] = {"objects":{}, "counts":{}, "next_id":0};));
        Regenerate();
    }

    void Regenerate() {
        pop = make_pop(r, pop_size, n_traits);
        eco_network = calc_competition(pop, eco_ea_fitness, settings);
        lex_network = calc_competition(pop, lexicase_fitness, settings);
        share_network = calc_competition(pop, sharing_fitness, settings);

        resources.resize(0);
        for (size_t i=0; i<n_traits; i++) {
            resources.push_back(emp::Resource(GetResourceInflow(), GetResourceInflow(), GetResourceOutflow()));
        }

    }

    void Update() {
        eco_network = calc_competition(pop, eco_ea_fitness, settings);
        lex_network = calc_competition(pop, lexicase_fitness, settings);
        share_network = calc_competition(pop, sharing_fitness, settings);        

        for (size_t i=0; i<n_traits; i++) {
            resources[i].Update();
        }

    }

    void SetPopSize(int s) {
        pop_size = s;
    }
    int GetPopSize() {
        return pop_size;
    }

    emp::vector<org_t> & GetPop() {
        return pop;
    }

    void SetNTraits(int n) {
        n_traits = n;
    }

    void SetSigmaShare(double s) {
        settings.SetSigmaShare(s);
    }
    double GetSigmaShare() {
        return SigmaShare::Get(settings);
    }

    void SetAlpha(double s) {
        settings.SetAlpha(s);
    }

    double GetAlpha() {
        return Alpha::Get(settings);
    }

    void SetResourceInflow(double s) {
        settings.SetResourceInflow(s);
    }

    double GetResourceInflow() {
        return ResourceInflow::Get(settings);
    }

    void SetResourceOutflow(double s) {
        settings.SetResourceOutflow(s);
    }

    double GetResourceOutflow() {
        return ResourceOutflow::Get(settings);
    }

    void SetCost(double s) {
        settings.SetCost(s);
    }

    double GetCost() {
        return Cost::Get(settings);
    }

    void SetCf(double s) {
        settings.SetCf(s);
    }

    double GetCf() {
        return Cf::Get(settings);
    }

    void SetNicheWidth(double s) {
        settings.SetNicheWidth(s);
    }

    double GetNicheWidth() {
        return NicheWidth::Get(settings);
    }

    void SetMaxScore(double s) {
        settings.SetMaxScore(s);
    }

    double GetMaxScore() {
        return MaxScore::Get(settings);
    }

    void SetMaxBonus(double s) {
        settings.SetMaxBonus(s);
    }

    double GetMaxBonus() {
        return MaxBonus::Get(settings);
    }

    void TournamentSelect() {
        emp::vector<org_t> new_pop;

        emp::vector<size_t> entries;
        for (size_t T = 0; T < pop_size; T++) {
            entries.resize(0);
            // Choose organisms for this tournament (with replacement!)
            for (size_t i=0; i < 2; i++) entries.push_back( r.GetUInt(pop_size) );

            double best_fit = emp::Sum(pop[entries[0]]);
            size_t best_id = entries[0];

            // Search for a higher fit org in the tournament.
            for (size_t i = 1; i < 2; i++) {
                const double cur_fit = emp::Sum(pop[entries[i]]);
                if (cur_fit > best_fit) {
                    best_fit = cur_fit;
                    best_id = entries[i];
                }
            }

            // Place the highest fitness into the next generation!
            new_pop.push_back(pop[best_id]);
        }
        pop = new_pop;
        Update();
    }

    void LexicaseSelect() {
        emp::vector<org_t> new_pop;

        std::map<org_t, int> genotype_counts;
        emp::vector<emp::vector<size_t>> genotype_lists;

        // Find all orgs with same genotype - we can dramatically reduce
        // fitness evaluations this way.
        for (size_t org_id = 0; org_id < pop_size; org_id++) {
            org_t gen = pop[org_id];
            if (emp::Has(genotype_counts, gen)) {
                genotype_lists[genotype_counts[gen]].push_back(org_id);
            } else {
                genotype_counts[gen] = genotype_lists.size();
                genotype_lists.emplace_back(emp::vector<size_t>{org_id});
            }
        }

        emp::vector<size_t> all_gens(genotype_lists.size()), cur_gens, next_gens;

        for (size_t i = 0; i < genotype_lists.size(); i++) {
            all_gens[i] = i;
        }

        emp::vector< emp::vector<double> > fitnesses(n_traits);

        for (size_t fit_id = 0; fit_id < n_traits; ++fit_id) {
            fitnesses[fit_id].resize(genotype_counts.size());
            int id = 0;
            for (auto & gen : genotype_lists) {
                fitnesses[fit_id][id] = gen[fit_id];
                id++;
            }
        }

        for (size_t repro = 0; repro < pop_size; ++repro) {
            emp::vector<size_t> order = emp::GetPermutation(r, n_traits);

            cur_gens = all_gens;  // Start with all of the organisms.
            int depth = -1;
            for (size_t fit_id : order) {
                depth++;

                double max_fit = fitnesses[fit_id][cur_gens[0]];
                next_gens.push_back(cur_gens[0]);
    
                for (size_t gen_id : cur_gens) {

                    const double cur_fit = fitnesses[fit_id][gen_id];
            
                    if (cur_fit > max_fit) {
                        max_fit = cur_fit;             // This is a the NEW maximum fitness for this function
                        next_gens.resize(0);           // Clear out orgs with former maximum fitness
                        next_gens.push_back(gen_id);   // Add this org as only one with new max fitness
                    }
                    else if (cur_fit == max_fit) {
                        next_gens.push_back(gen_id);   // Same as cur max fitness -- save this org too.
                    }
                }
                // Make next_orgs into new cur_orgs; make cur_orgs allocated space for next_orgs.
                std::swap(cur_gens, next_gens);
                next_gens.resize(0);

                if (cur_gens.size() == 1) break;  // Stop if we're down to just one organism.
            }

            // Place a random survivor (all equal) into the next generation!
            emp_assert(cur_gens.size() > 0, cur_gens.size(), n_traits, all_gens.size());
            size_t options = 0;
            for (size_t gen : cur_gens) {
                options += genotype_lists[gen].size();
            }
            size_t winner = r.GetUInt(options);
            int repro_id = -1;

            for (size_t gen : cur_gens) {
                if (winner < genotype_lists[gen].size()) {
                    repro_id = (int) genotype_lists[gen][winner];
                    break;
                }
                winner -= genotype_lists[gen].size();
            }
            emp_assert(repro_id != -1, repro_id, winner, options);

            emp::vector<size_t> used = emp::Slice(order, 0, depth+1);
            new_pop.push_back(pop[repro_id]);
        }

        pop = new_pop;
        Update();
    }

    void ResourceSelect() {
        emp::vector<org_t> new_pop;

       emp::vector<double> base_fitness(pop_size);
       emp::vector< emp::vector<double> > extra_fitnesses(n_traits);
       for (size_t i=0; i < n_traits; i++) {
         extra_fitnesses[i].resize(pop_size);
       }

       for (size_t org_id = 0; org_id < pop_size; org_id++) {
            base_fitness[org_id] = 0;

            for (size_t ex_id = 0; ex_id < n_traits; ex_id++) {

                resources[ex_id].Inc(resources[ex_id].GetInflow()/pop_size);
                double cur_fit = pop[org_id][ex_id]/GetMaxScore();

                if (cur_fit > GetNicheWidth()) {
                    cur_fit = emp::Pow(cur_fit, 2.0);
                    cur_fit *= GetCf()*(resources[ex_id].GetAmount()-GetCost());
                    cur_fit -= GetCost();
                } else {
                    cur_fit = 0;
                }

                cur_fit = std::min(cur_fit, GetMaxBonus());
                extra_fitnesses[ex_id][org_id] = emp::Pow2(cur_fit);
                base_fitness[org_id] *= emp::Pow2(cur_fit);
                resources[ex_id].Dec(std::abs(cur_fit));
            }
       }

       emp::vector<size_t> entries;
       for (size_t T = 0; T < pop_size; T++) {
         entries.resize(0);

         for (size_t i=0; i<2; i++) entries.push_back( r.GetUInt(pop_size) ); // Allows replacement!

        double best_fit = base_fitness[entries[0]];
        size_t best_id = entries[0];

        // Search for a higher fit org in the tournament.
        for (size_t i = 1; i < 2; i++) {
            const double cur_fit = base_fitness[entries[i]];
            if (cur_fit > best_fit) {
                best_fit = cur_fit;
                best_id = entries[i];
            }
        }

         // Place the highest fitness into the next generation!
         new_pop.push_back(pop[best_id]);
       }


        for (size_t i=0; i<n_traits; i++) {
            resources[i].Dec(GetResourceInflow()); // Inflow is about to get added again
        }
        pop = new_pop;
        Update();

    }

    void SharingSelect() {
        emp::vector<double> fits;
        emp::vector<org_t> new_pop;

        double sharing_threshold = GetSigmaShare();
        double alpha = GetAlpha();

        for (org_t & org1 : pop) {
            double niche_count = 0;
            for (org_t & org2 : pop) {
                double dist = emp::EuclideanDistance(org1, org2);
                niche_count += std::max(1.0 - std::pow(dist/sharing_threshold, alpha), 0.0);            
            }
            fits.push_back(emp::Sum(org1)/niche_count);
        }

        emp::vector<size_t> entries;
        for (size_t T = 0; T < pop_size; T++) {
            entries.resize(0);
            // Choose organisms for this tournament (with replacement!)
            for (size_t i=0; i < 2; i++) entries.push_back( r.GetUInt(pop_size) );

            double best_fit = fits[entries[0]];
            size_t best_id = entries[0];

            // Search for a higher fit org in the tournament.
            for (size_t i = 1; i < 2; i++) {
                const double cur_fit = fits[entries[i]];
                if (cur_fit > best_fit) {
                    best_fit = cur_fit;
                    best_id = entries[i];
                }
            }

            // Place the highest fitness into the next generation!
            new_pop.push_back(pop[best_id]);
        }
        pop = new_pop;
        Update();
    }

};


#endif