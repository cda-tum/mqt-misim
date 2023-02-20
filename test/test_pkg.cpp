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
  auto dd = std::make_unique<dd::MDDPackage>(2, std::vector<std::size_t>{3, 2});
  EXPECT_EQ(dd->qregisters(), 2);

  auto hGate0 = dd->makeGateDD<dd::TritMatrix>(dd::H3, 2, 0);

  ASSERT_TRUE(dd->getValueByPath(hGate0, "00")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "10")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "20")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "30")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "40")
                  .approximatelyEquals(dd::ComplexValue{-0.28867513, 0.5}));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "50")
                  .approximatelyEquals(dd::ComplexValue{-0.28867513, -0.5}));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "60")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "70")
                  .approximatelyEquals(dd::ComplexValue{-0.28867513, -0.5}));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "80")
                  .approximatelyEquals(dd::ComplexValue{-0.28867513, 0.5}));

  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "01").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "11").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "21").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "31").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "41").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "51").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "61").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "71").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "81").approximatelyEquals(dd::COMPLEX_ZERO));

  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "02").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "12").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "22").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "32").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "42").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "52").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "62").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "72").approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(
      dd->getValueByPath(hGate0, "82").approximatelyEquals(dd::COMPLEX_ZERO));

  ASSERT_TRUE(dd->getValueByPath(hGate0, "03")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "13")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "23")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "33")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "43")
                  .approximatelyEquals(dd::ComplexValue{-0.28867513, 0.5}));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "53")
                  .approximatelyEquals(dd::ComplexValue{-0.28867513, -0.5}));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "63")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "73")
                  .approximatelyEquals(dd::ComplexValue{-0.28867513, -0.5}));
  ASSERT_TRUE(dd->getValueByPath(hGate0, "83")
                  .approximatelyEquals(dd::ComplexValue{-0.28867513, 0.5}));

  // case with higher target
  auto hGate1 = dd->makeGateDD<dd::GateMatrix>(dd::H2, 2, 1);
  ASSERT_EQ(dd->getValueByPath(hGate1, "00"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "10"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "20"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "30"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "40"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "50"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "60"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "70"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "80"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));

  ASSERT_EQ(dd->getValueByPath(hGate1, "01"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "11"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "21"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "31"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "41"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "51"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "61"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "71"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "81"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));

  ASSERT_EQ(dd->getValueByPath(hGate1, "02"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "12"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "22"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "32"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "42"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "52"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "62"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "72"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "82"),
            (dd::ComplexValue{dd::SQRT2_2, 0}));

  ASSERT_EQ(dd->getValueByPath(hGate1, "03"),
            (dd::ComplexValue{-dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "13"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "23"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "33"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "43"),
            (dd::ComplexValue{-dd::SQRT2_2, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "53"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "63"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "73"), (dd::ComplexValue{0, 0}));
  ASSERT_EQ(dd->getValueByPath(hGate1, "83"),
            (dd::ComplexValue{-dd::SQRT2_2, 0}));

  dd::Controls controls1{{1, 1}};
  auto controlledH3 = dd->makeGateDD<dd::TritMatrix>(dd::H3, 2, controls1, 0);

  ASSERT_TRUE(dd->getValueByPath(controlledH3, "00")
                  .approximatelyEquals(dd::COMPLEX_ONE));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "10")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "20")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "30")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "40")
                  .approximatelyEquals(dd::COMPLEX_ONE));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "50")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "60")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "70")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "80")
                  .approximatelyEquals(dd::COMPLEX_ONE));

  ASSERT_TRUE(dd->getValueByPath(controlledH3, "01")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "11")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "21")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "31")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "41")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "51")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "61")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "71")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "81")
                  .approximatelyEquals(dd::COMPLEX_ZERO));

  ASSERT_TRUE(dd->getValueByPath(controlledH3, "02")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "12")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "22")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "32")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "42")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "52")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "62")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "72")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "82")
                  .approximatelyEquals(dd::COMPLEX_ZERO));

  ASSERT_TRUE(dd->getValueByPath(controlledH3, "03")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "13")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "23")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "33")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "43")
                  .approximatelyEquals(dd::ComplexValue{-0.28867513, 0.5}));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "53")
                  .approximatelyEquals(dd::ComplexValue{-0.28867513, -0.5}));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "63")
                  .approximatelyEquals(dd::COMPLEX_SQRT3_3));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "73")
                  .approximatelyEquals(dd::ComplexValue{-0.28867513, -0.5}));
  ASSERT_TRUE(dd->getValueByPath(controlledH3, "83")
                  .approximatelyEquals(dd::ComplexValue{-0.28867513, 0.5}));

  dd::Controls controls0{{0, 1}};
  auto controlled0H2 = dd->makeGateDD<dd::GateMatrix>(dd::H2, 2, controls0, 1);

  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "00")
                  .approximatelyEquals(dd::COMPLEX_ONE));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "10")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "20")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "30")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "40")
                  .approximatelyEquals(dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "50")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "60")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "70")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "80")
                  .approximatelyEquals(dd::COMPLEX_ONE));

  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "01")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "11")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "21")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "31")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "41")
                  .approximatelyEquals(dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "51")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "61")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "71")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "81")
                  .approximatelyEquals(dd::COMPLEX_ZERO));

  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "02")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "12")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "22")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "32")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "42")
                  .approximatelyEquals(dd::ComplexValue{dd::SQRT2_2, 0}));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "52")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "62")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "72")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "82")
                  .approximatelyEquals(dd::COMPLEX_ZERO));

  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "03")
                  .approximatelyEquals(dd::COMPLEX_ONE));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "13")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "23")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "33")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "43")
                  .approximatelyEquals(dd::ComplexValue{-dd::SQRT2_2, 0}));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "53")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "63")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "73")
                  .approximatelyEquals(dd::COMPLEX_ZERO));
  ASSERT_TRUE(dd->getValueByPath(controlled0H2, "83")
                  .approximatelyEquals(dd::COMPLEX_ONE));

  using namespace dd::literals;
}

TEST(DDPackageTest, Identity) {
  auto dd =
      std::make_unique<dd::MDDPackage>(4, std::vector<std::size_t>{2, 3, 2, 4});
  EXPECT_EQ(dd->qregisters(), 4);
  EXPECT_TRUE(dd->makeIdent(0).isOneTerminal());
  EXPECT_TRUE(dd->makeIdent(0, -1).isOneTerminal());

  auto id3 = dd->makeIdent(3);
  EXPECT_EQ(dd->makeIdent(0, 2), id3);
  const auto& table = dd->getIdentityTable();
  EXPECT_NE(table[2].nextNode, nullptr);

  auto id2 = dd->makeIdent(0, 1);  // should be found in idTable
  EXPECT_EQ(dd->makeIdent(2), id2);

  auto id4 = dd->makeIdent(0, 3);  // should use id3 and extend it
  EXPECT_EQ(dd->makeIdent(0, 3), id4);
  EXPECT_NE(table[3].nextNode, nullptr);

  auto idCached = dd->makeIdent(4);
  EXPECT_EQ(id4, idCached);
}

TEST(DDPackageTest, Multiplication) {
  auto dd = std::make_unique<dd::MDDPackage>(2, std::vector<std::size_t>{2});
  EXPECT_EQ(dd->qregisters(), 2);

  auto hGate1 = dd->makeGateDD<dd::GateMatrix>(dd::H2, 1, 0);

  auto zeroState = dd->makeZeroState(1);
  auto hState = dd->multiply(hGate1, zeroState);
  auto x = 0;
  /*
  ASSERT_EQ(dd->fidelity(zeroState, oneState), 0.0);
  // repeat the same calculation - triggering compute table hit
  ASSERT_EQ(dd->fidelity(zeroState, oneState), 0.0);
  ASSERT_NEAR(dd->fidelity(zeroState, hState), 0.5,
  dd::ComplexTable<>::tolerance()); ASSERT_NEAR(dd->fidelity(oneState, hState),
  0.5, dd::ComplexTable<>::tolerance());
   */
}