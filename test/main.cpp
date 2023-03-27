#include "dd/MDDPackage.hpp"

#include <memory>
#include <vector>

int main() { // NOLINT(bugprone-exception-escape)
    auto mdd = std::make_unique<dd::MDDPackage>(1, std::vector<std::size_t>{5});

    auto rxy = mdd->makeGateDD<dd::QuintMatrix>(dd::RXY5(dd::PI_4, dd::PI_4, 0, 1), 1, 0);

    // Gates
    //auto h3Gate = mdd->makeGateDD<dd::QuintMatrix>(dd::H5(), 1, 0);
    mdd->getVectorizedMatrix(rxy);

    auto oneState = mdd->makeBasisState(1, {1});

    mdd->printVector(oneState);

    return 0;
}
