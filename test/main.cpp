#include <iostream>

#include "dd/MDDPackage.hpp"

int main() {
  auto mdd = dd::MDDPackage(3, {2, 2, 2});
  auto zeroState = mdd.makeZeroState(3);  // number of particles in the rack
  mdd.printVector(zeroState);
  auto id = mdd.makeIdent(3);

  auto x = mdd.getVector(zeroState);
  auto hGate = mdd.makeGateDD<dd::GateMatrix>(dd::Hmat, 2, 1);

  return 0;
}
