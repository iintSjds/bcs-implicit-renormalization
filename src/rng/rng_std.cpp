#include "rng.hpp"

// using standard c++ library
#include <random>

namespace RNG{
    std::random_device dev;
    std::mt19937_64 generator(dev());
    std::uniform_real_distribution<double> dist(0,1);
    double surn(unsigned long seed){
	generator.seed(seed);
	return dist(generator);
    }    
    double urn(){
	return dist(generator);
    }
}
