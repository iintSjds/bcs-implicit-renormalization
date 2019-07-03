#include <iostream>
#include "rng.hpp"

int main(){
    for(int i=0;i<5;i++)
	std::cout<<RNG::urn()<<"\t";
    std::cout<<std::endl;
    for(int i=0;i<5;i++)
	std::cout<<RNG::surn(1)<<"\t";
    std::cout<<std::endl;

    return 0;
}
