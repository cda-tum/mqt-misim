#include "dd/MDDPackage.hpp"

#include <memory>
#include <vector>

int main() { // NOLINT(bugprone-exception-escape)
    auto mdd = std::make_unique<dd::MDDPackage>(3, std::vector<std::size_t>{2, 2, 3});

    auto oneState = mdd->makeBasisState(3, {1, 1, 1});

    mdd->printVector(oneState);

    return 0;
}
