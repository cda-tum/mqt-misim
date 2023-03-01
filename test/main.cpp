#include "dd/MDDPackage.hpp"

#include <memory>
#include <vector>

int main() { // NOLINT(bugprone-exception-escape)
    auto mdd =
            std::make_unique<dd::MDDPackage>(3, std::vector<std::size_t>{2, 2, 3});
    // auto zeroState = mdd.makeZeroState(3);  // number of particles in the rack
    auto oneState = mdd->makeBasisState(3, {1, 1, 2});
    mdd->printVector(oneState);
    // auto id = mdd.makeIdent(3);

    // auto x = mdd.getVector(zeroState);
    // auto U3 = dd::U3mat(dd::PI_4, -dd::PI_2, dd::PI_2);

    return 0;
}
