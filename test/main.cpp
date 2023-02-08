#include <iostream>

#include "dd/MDDPackage.hpp"

int main() {
  auto mdd = dd::MDDPackage(2, {2, 2});
  auto zeroState = mdd.makeZeroState(2);  // number of particles in the rack
  mdd.printVector(zeroState);
  auto id = mdd.makeIdent(2);

  auto x = mdd.getVector(zeroState);
  auto hGate = mdd.makeGateDD<dd::GateMatrix>(dd::Hmat, 1, 0);

  return 0;
}
