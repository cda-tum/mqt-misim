#include <iostream>

#include "dd/MDDPackage.hpp"

int main() {
  auto mdd = dd::MDDPackage(3, {2, 2, 3});
  // auto zeroState = mdd.makeZeroState(3);  // number of particles in the rack
  auto oneState = mdd.makeBasisState(3, {1, 1, 2});
  mdd.printVector(oneState);
  // auto id = mdd.makeIdent(3);

  // auto x = mdd.getVector(zeroState);
  // auto U3 = dd::U3mat(dd::PI_4, -dd::PI_2, dd::PI_2);

  return 0;
}
