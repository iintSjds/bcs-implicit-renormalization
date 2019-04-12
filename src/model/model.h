#ifndef _SIMPLE_SC_MODEL_
#define _SIMPLE_SC_MODEL_

#include <vector>

/*  ___  ___   __  __  ___  ___  ___ _    ___  */
/* / __|/ __| |  \/  |/ _ \|   \| __| |  / __| */
/* \__ \ (__  | |\/| | (_) | |) | _|| |__\__ \ */
/* |___/\___| |_|  |_|\___/|___/|___|____|___/ */

// currently define simple models

// Model should include all information specified by the model,
// including coordinates, meshgrid, functions(Gamma and GG), integration methods, and inner product, etc.


class Model{
 public:
    int N1,N2; // total number of points of inner/outer region;   N1+N2
    std::vector<double> delta1,delta2,freq1,freq2;
    //TBA
    
    
};
