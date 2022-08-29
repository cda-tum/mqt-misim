#include <iostream>
#include "include/dd/MDDPackage.hpp"

int main() {
    //dd::MDDPackage mdd();
    auto mdd = dd::MDDPackage( 2, {2,3} );
    auto zero_state = mdd.makeZeroState(1);
    mdd.printVector(zero_state);
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
