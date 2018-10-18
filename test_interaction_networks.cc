#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include "interaction_networks.h"
// Controller t;

TEST_CASE("MakePop", "[helpers]") {
    emp::Random r;
    emp::vector<org_t> pop = make_pop(r, 20, 10);
    REQUIRE(pop.size() == 20);
    CHECK(pop[0].size() == 10);
}

TEST_CASE("FindHighest", "[helpers]") {
    emp::vector<org_t> v({{1,2,3}, {2, 1, 3}, {1,3,2}});
    auto highest = FindHighest(v, 0);
    CHECK(highest.size() == 1);
    emp::vector<int> res = highest[0];

    for (size_t i = 0; i < res.size(); i++) {
        CHECK( res[i] == v[1][i]);
    }

    highest = FindHighest(v, 1);
    CHECK(highest.size() == 1);
    res = highest[0];

    for (size_t i = 0; i < res.size(); i++) {
        CHECK( res[i] == v[2][i]);
    }

    highest = FindHighest(v, 2);
    CHECK(highest.size() == 2);
}

TEST_CASE("Lexicase", "[selection_schemes]") {
    emp::vector<org_t> pop({{3,0,0}, {0, 3, 0}, {0,0,3}});
    fit_map_t fits = lexicase_fitness(pop, DEFAULT);
    for (auto & o : fits) {
        CHECK(o.second == Approx(.33333333333));
    }

    pop = emp::vector<org_t>({{3,3,3}, {0, 1, 2}, {2,1,1}});
    fits = lexicase_fitness(pop, DEFAULT);
    CHECK(fits[{3,3,3}] == 1);
    CHECK(fits[{0,1,2}] == 0);
    CHECK(fits[{2,1,1}] == 0);

    pop = emp::vector<org_t>({{3,3,3}, {3, 1, 2}, {2,1,1}});
    fits = lexicase_fitness(pop, DEFAULT);
    CHECK(fits[{3,3,3}] == 1);
    CHECK(fits[{3,1,2}] == 0);
    CHECK(fits[{2,1,1}] == 0);

    pop = emp::vector<org_t>({{3,3,3}, {3, 3, 3}, {2,1,1}});
    fits = lexicase_fitness(pop, DEFAULT);
    CHECK(fits[{3,3,3}] == Approx(.5));
    CHECK(fits[{2,1,1}] == 0);

    pop = emp::vector<org_t>({{3,1,2}, {1, 1, 2}, {2,1,1}});
    fits = lexicase_fitness(pop, DEFAULT);
    CHECK(fits[{3,1,2}] == 1);
    CHECK(fits[{1,1,2}] == 0);
    CHECK(fits[{2,1,1}] == 0);

    pop = emp::vector<org_t>({{3,1,2}, {1, 3, 2}, {2,3,1}});
    fits = lexicase_fitness(pop, DEFAULT);
    CHECK(fits[{3,1,2}] == Approx(.5));
    CHECK(fits[{1,3,2}] == Approx(.333333));
    CHECK(fits[{2,3,1}] == Approx(.16666667));

    pop = emp::vector<org_t>({{3,1,2,1,1}, {1, 3, 2,1,1}, {2,3,1,1,1}});
    fits = lexicase_fitness(pop, DEFAULT);
    CHECK(fits[{3,1,2,1,1}] == Approx(.5));
    CHECK(fits[{1,3,2,1,1}] == Approx(.333333));
    CHECK(fits[{2,3,1,1,1}] == Approx(.16666667));


}

TEST_CASE("Fitness sharing", "[selection_schemes]") {
    emp::vector<org_t> pop = emp::vector<org_t>({{3,3,3}, {3, 3, 3}, {3,3,3}});
    all_attrs settings = DEFAULT;
    fit_map_t fits = sharing_fitness(pop, settings);
    CHECK(fits[{3,3,3}] == Approx(.3333333));

    pop = emp::vector<org_t>({{3,1,2,1,1}, {1, 3, 2,1,1}, {2,3,1,1,1}});
    fits = sharing_fitness(pop, DEFAULT);
    CHECK(fits[{3,1,2,1,1}] == Approx(.666666));
    CHECK(fits[{1,3,2,1,1}] == Approx(.333333));
    CHECK(fits[{2,3,1,1,1}] == Approx(0));

    pop = emp::vector<org_t>({{10,1}, {1, 10}, {1,1}});
    fits = sharing_fitness(pop, DEFAULT);
    CHECK(fits[{10,1}] == Approx(.5));
    CHECK(fits[{1,10}] == Approx(.5));
    CHECK(fits[{1,1}] == Approx(0));

}

TEST_CASE("Eco-EA", "[selection_schemes]") {
    emp::vector<org_t> pop = emp::vector<org_t>({{3,3,3}, {3, 3, 3}, {3,3,3}});
    all_attrs settings = DEFAULT;
    fit_map_t fits = eco_ea_fitness(pop, settings);
    CHECK(fits[{3,3,3}] == Approx(.3333333));

    pop = emp::vector<org_t>({{3,1,2,1,1}, {1, 3, 2,1,1}, {2,3,1,1,1}});
    fits = eco_ea_fitness(pop, DEFAULT);
    CHECK(fits[{3,1,2,1,1}] == Approx(.666666));
    CHECK(fits[{1,3,2,1,1}] == Approx(.1666666));
    CHECK(fits[{2,3,1,1,1}] == Approx(.16666667));

    pop = emp::vector<org_t>({{10,1,2,1,1}, {1, 3, 2,1,1}, {2,3,1,1,1}, {2,1,1,1,1}});
    fits = eco_ea_fitness(pop, DEFAULT);
    CHECK(fits[{10,1,2,1,1}] == Approx(.5));
    CHECK(fits[{1,3,2,1,1}] == Approx(0.0833335));
    CHECK(fits[{2,3,1,1,1}] == Approx(0.0833335));
    CHECK(fits[{2,1,1,1,1}] == Approx(.3333333));

}

TEST_CASE("Calc competition", "[helpers]") {
    emp::vector<org_t> pop = emp::vector<org_t>({{1,3,1}, {3, 1, 1}, {1,1,3}});
    // all_attrs settings = DEFAULT;
    std::function<fit_map_t(emp::vector<org_t>&, all_attrs)> test_fun = [](emp::vector<org_t> & pop, all_attrs attrs=DEFAULT) {
        fit_map_t base_fit_map;
        for (org_t & org : pop) {
            base_fit_map[org] = 1.0;
        }
        return base_fit_map;
    };

    emp::WeightedGraph g = calc_competition(pop, test_fun);
    auto weights = g.GetWeights();
    for (auto vec : weights) {
        for (auto val : vec) {
            CHECK(val == 0); 
        }
    }

    test_fun = [](emp::vector<org_t> & pop, all_attrs attrs=DEFAULT) {
        fit_map_t base_fit_map;
        for (org_t & org : pop) {
            base_fit_map[org] = 1.0;
        }

        if (emp::Has(pop, {1,3,1})) {
            base_fit_map[{3,1,1}] = 0;
        }
        return base_fit_map;
    };

    g = calc_competition(pop, test_fun);
    CHECK(g.GetWeight(0,1) == -1);
}

// TEST_CASE("Controller", "[controller]") {
//     Controller c;
//     c.SetPopSize(50);
//     CHECK(c.GetPopSize() == 50);
//     CHECK(c.GetPop().size() == 10);
//     c.Regenerate();
//     CHECK(c.GetPop().size() == 50);

//     c.SetNTraits(9);
//     c.SetSigmaShare(3);
//     c.SetAlpha(.1);
//     c.SetCost(2);
//     c.SetCf(.01);
//     c.SetNicheWidth(4);
//     c.SetMaxScore(100);

//     c.Regenerate();

//     CHECK(c.GetPop()[0].size() == 9);
//     CHECK(c.GetSigmaShare() == 3);
//     CHECK(Alpha::Get(c.settings) == .1);
//     CHECK(Cost::Get(c.settings) == 2);
//     CHECK(Cf::Get(c.settings) == .01);
//     CHECK(NicheWidth::Get(c.settings) == 4);
//     CHECK(MaxScore::Get(c.settings) == 100);
// }