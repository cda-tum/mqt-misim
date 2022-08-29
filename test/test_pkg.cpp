//
// Created by k3vn on 27.08.22.
//
//#include "dd/Export.hpp"
/*
#include "include/dd/GateMatrixDefinitions.hpp"
#include "include/dd/MDDPackage.hpp"

//#include "gtest/gtest.h"
#include <iomanip>
#include <memory>
#include <random>
#include <sstream>

using namespace dd::literals;

TEST(DDPackageTest, TrivialTest) {

auto dd = std::make_unique<dd::Package>(2);

EXPECT_EQ(dd->qubits(), 2);

auto x_gate = dd->makeGateDD(dd::Xmat, 1, 0);


ASSERT_EQ(dd->getValueByPath(h_gate, "0"), (dd::ComplexValue{dd::SQRT2_2, 0}));

auto zero_state = dd->makeZeroState(1);
auto h_state    = dd->multiply(h_gate, zero_state);
auto one_state  = dd->multiply(x_gate, zero_state);

ASSERT_EQ(dd->fidelity(zero_state, one_state), 0.0);
// repeat the same calculation - triggering compute table hit
ASSERT_EQ(dd->fidelity(zero_state, one_state), 0.0);
ASSERT_NEAR(dd->fidelity(zero_state, h_state), 0.5, dd::ComplexTable<>::tolerance());
ASSERT_NEAR(dd->fidelity(one_state, h_state), 0.5, dd::ComplexTable<>::tolerance());
}
*/