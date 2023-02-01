#include "dd/MDDPackage.hpp"

#include <iostream>

int main() {
    //dd::MDDPackage mdd();
    auto mdd = dd::MDDPackage(3, {2, 3, 3});
    auto zero_state = mdd.makeZeroState(3); // number of particles in the rack
    mdd.printVector(zero_state);
    auto id = mdd.makeIdent(3);

    auto x = mdd.getVector(zero_state);

    return 0;
}
