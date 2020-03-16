#include <vector>
#include <algorithm>
#include <random>

std::vector<int> make_population(int nbr_pts)
{
    std::vector<int> pop(nbr_pts);
    std::iota(std::begin(pop), std::end(pop), 0);
    return pop;
}

std::vector<int> randperm(int N, int nbr_pts)
{
    static std::vector<int> population = make_population(nbr_pts);
    static std::mt19937 rng(std::random_device{}());

    // return the first N numbers of the randomly shuffled vector
    std::shuffle(std::begin(population), std::end(population), rng);
    return {std::begin(population), std::begin(population)+N};
}
