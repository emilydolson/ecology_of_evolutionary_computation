#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include "interaction_networks.h"

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

