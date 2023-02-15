#include <iostream>

#include "dd/MDDPackage.hpp"

int main() {
  auto mdd = dd::MDDPackage(1, {2});
  // auto zeroState = mdd.makeZeroState(3);  // number of particles in the rack

  // mdd.printVector(zeroState);
  // auto id = mdd.makeIdent(3);

  // auto x = mdd.getVector(zeroState);
  auto U3 = dd::U3mat(dd::PI_4, -dd::PI_2, dd::PI_2);

  return 0;
}
