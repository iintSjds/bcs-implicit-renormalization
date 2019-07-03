#include "rng.hpp"

//using trng yarn2 generator
#include <trng/yarn2.hpp>
#include <trng/uniform01_dist.hpp>

namespace RNG{

    trng::yarn2 generator;
    trng::uniform01_dist<> dist;
    double surn(unsigned long seed){
	generator.seed(seed);
	return dist(generator);
    }
    double urn(){
	return dist(generator);
    }
}
