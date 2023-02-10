#include <cstddef>

#include "dd/MDDPackage.hpp"
#include "gtest/gtest.h"

using namespace dd::literals;

TEST(DDPackageTest, RequestInvalidPackageSize) {
  EXPECT_THROW(auto dd = std::make_unique<dd::MDDPackage>(
                   std::numeric_limits<dd::QuantumRegister>::max() + 2,
                   std::vector<std::size_t>(
                       2, std::numeric_limits<dd::QuantumRegister>::max() + 2)),
               std::invalid_argument);
}

TEST(DDPackageTest, TrivialTest) {
  auto dd = std::make_unique<dd::MDDPackage>(2, std::vector<std::size_t>{2, 2});
  EXPECT_EQ(dd->qregisters(), 2);

  // case with higher target
  auto hGate1 = dd->makeGateDD<dd::GateMatrix>(dd::Hmat, 2, 1);
  ASSERT_EQ(dd->getValueByPath(hGate1, "00"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "01"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "02"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "03"),
            (dd::ComplexValue{-dd::SQRT2_2, 0}));

  ASSERT_EQ(dd->getValueByPath(hGate1, "10"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "11"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "12"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "13"), (dd::ComplexValue{0, 0}));

  ASSERT_EQ(dd->getValueByPath(hGate1, "20"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "21"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "22"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "23"), (dd::ComplexValue{0, 0}));

  ASSERT_EQ(dd->getValueByPath(hGate1, "30"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "31"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "32"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "33"),
            (dd::ComplexValue{-dd::SQRT2_2, 0}));

  // case with lower target
  auto hGate0 = dd->makeGateDD<dd::GateMatrix>(dd::Hmat, 2, 0);
  ASSERT_EQ(dd->getValueByPath(hGate0, "00"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate0, "01"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate0, "02"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate0, "03"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));

  ASSERT_EQ(dd->getValueByPath(hGate0, "10"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate0, "11"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate0, "12"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate0, "13"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));

  ASSERT_EQ(dd->getValueByPath(hGate0, "20"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate0, "21"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate0, "22"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate0, "23"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));

  ASSERT_EQ(dd->getValueByPath(hGate0, "30"),
            (dd::ComplexValue{-dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate0, "31"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate0, "32"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate0, "33"),
            (dd::ComplexValue{-dd::SQRT2_2, 0}));

  using namespace dd::literals;
}
/*
TEST(DDPackageTest, Identity) {
  auto dd = std::make_unique<dd::MDDPackage>(4,
std::vector<std::size_t>{2,2,2,2});

  EXPECT_TRUE(dd->makeIdent(0).isOneTerminal());
  EXPECT_TRUE(dd->makeIdent(0, -1).isOneTerminal());

  auto id3 = dd->makeIdent(3);
  EXPECT_EQ(dd->makeIdent(0, 2), id3);
  const auto& table = dd->getIdentityTable();
  EXPECT_NE(table[2].nextNode, nullptr);

  auto id2 = dd->makeIdent(0, 1); // should be found in idTable
  EXPECT_EQ(dd->makeIdent(2), id2);

  auto id4 = dd->makeIdent(0, 3); // should use id3 and extend it
  EXPECT_EQ(dd->makeIdent(0, 3), id4);
  EXPECT_NE(table[3].nextNode, nullptr);

  auto idCached = dd->makeIdent(4);
  EXPECT_EQ(id4, idCached);
}
 */