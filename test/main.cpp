#include "dd/MDDPackage.hpp"

#include <memory>
#include <vector>

int main() { // NOLINT(bugprone-exception-escape)
    auto mdd = std::make_unique<dd::MDDPackage>(2, std::vector<std::size_t>{5, 5});

    //CSUM(QuantumRegisterCount n, QuantumRegister cReg, QuantumRegister target)
    auto h5Gate = mdd->makeGateDD<dd::QuintMatrix>(dd::H5(), 2, 1);
    //auto csum = mdd->CSUM(2, 1, 0);

    auto mul = mdd->makeGateDD<dd::QuintMatrix>(dd::I5, 2, 1);
    mul      = mdd->multiply(mul, h5Gate);
    mul      = mdd->multiply(mul, h5Gate);
    mul      = mdd->multiply(mul, h5Gate);
    mul      = mdd->multiply(mul, h5Gate);

    //mul = mdd->multiply(mul, h5Gate);
    mdd->getVectorizedMatrix(mul);

    //auto oneState = mdd->makeBasisState(1, {1});

    //mdd->printVector(oneState);

    return 0;
}
