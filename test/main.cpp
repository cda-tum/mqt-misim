#include <iostream>

#include "dd/MDDPackage.hpp"

int main() {
  auto mdd = dd::MDDPackage(3, {2, 2, 3});

    auto oneState = mdd.makeBasisState(3, {1, 1, 2});

  mdd.printVector(oneState);

  return 0;
}
