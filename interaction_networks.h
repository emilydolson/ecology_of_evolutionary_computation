#ifndef _INTERACTION_NETWORKS_H_
#define _INTERACTION_NETWORKS_H_

#include <functional>
#include <map>
#include <limits>

#include "web/init.h"

#include "base/vector.h"
#include "base/assert.h"
#include "tools/Random.h"
#include "tools/vector_utils.h"
#include "tools/string_utils.h"
#include "tools/math.h"
#include "tools/distances.h"
#include "tools/attrs.h"
#include "tools/Graph.h"

#include <permutations.hpp>

using org_t = emp::vector<int>;
using fit_map_t = std::map<org_t, double>;


DEFINE_ATTR(SigmaShare);
DEFINE_ATTR(Alpha);
DEFINE_ATTR(Cost);
DEFINE_ATTR(Cf);
DEFINE_ATTR(NicheWidth);
DEFINE_ATTR(MaxScore);

constexpr auto DEFAULT{MakeAttrs(SigmaShare(8.0),
                                 Alpha(1.0),
                                 Cost(1.0),
                                 Cf(.0025),                                 
                                 NicheWidth(3.0),
                                 MaxScore(10.0))};

using all_attrs = emp::tools::Attrs<typename SigmaShare::value_t<double>, typename Alpha::value_t<double>, 
                        typename Cost::value_t<double>, typename Cf::value_t<double>,
                        typename NicheWidth::value_t<double>, typename MaxScore::value_t<double>>;


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

void TraverseDecisionTree(fit_map_t & fit_map, emp::vector<org_t> & pop, emp::vector<int> axes) {
    if (axes.size() == 1) {
        emp::vector<org_t> winners = FindHighest(pop, axes[0]);
        for (org_t & org : winners) {
            fit_map[org]+=1.0/(double)winners.size();
        }
        return;
    }

    for (int ax : axes) {
        emp::vector<org_t> winners = FindHighest(pop, ax);
        if (winners.size() == 1) { // Not a tie
            fit_map[winners[0]] += (double)emp::Factorial(axes.size() - 1);
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
        fit_map[org] /= (double)emp::Factorial(n_funs); // convert to proportion of "islands"

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
        double res = 2000;
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
                base_fit_map[org] *= emp::Pow2(Cf::Get(attrs)*res*emp::Pow(org[axis]/MaxScore::Get(attrs),2.0) - Cost::Get(attrs));
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
        emp::vector<org_t> curr = pop;
        for (int & ax : curr[i]) {
            ax = 0; // Replace org with null org so pop size stays same
        }

        fit_map_t new_fits = fit_fun(curr, attrs);
        for (size_t j = 0; j < pop.size(); j++ ) {
            if (i == j) {continue;}
            // std::cout << fitnesses[curr[j]] << " " << new_fits[curr[j]] << " " << fitnesses[curr[j]] - new_fits[curr[j]]<< std::endl;
            double effect = fitnesses[curr[j]] - new_fits[curr[j]];
            effects.AddEdge(i, j, effect);
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

    void SetCost(double s) {
        settings.SetCost(s);
    }

    void SetCf(double s) {
        settings.SetCf(s);
    }

    void SetNicheWidth(double s) {
        settings.SetNicheWidth(s);
    }
    
    void SetMaxScore(double s) {
        settings.SetMaxScore(s);
    }


};


#endif