#include <iostream>
#include "include/dd/MDDPackage.hpp"

int main() {
    //dd::MDDPackage mdd();
    auto mdd = dd::MDDPackage( 2, {2,3} );
    auto zero_state = mdd.makeZeroState(2); // number of particles in the rack
    mdd.printVector(zero_state);

    auto x = mdd.getVector(zero_state);

    return 0;
}
