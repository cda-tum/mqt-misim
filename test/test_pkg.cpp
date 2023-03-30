#include "dd/MDDPackage.hpp"

#include "gtest/gtest.h"
#include <cstddef>

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

    dd::ComplexValue h3plus  = dd::COMPLEX_SQRT3_3 * dd::ComplexValue{std::cos(2. * dd::PI / 3.), std::sin(2. * dd::PI / 3.)};
    dd::ComplexValue h3minus = dd::COMPLEX_SQRT3_3 * dd::ComplexValue{std::cos(4. * dd::PI / 3.), std::sin(4. * dd::PI / 3.)};

    auto hGate0 = dd->makeGateDD<dd::TritMatrix>(dd::H3(), 2, 0);

    ASSERT_TRUE(dd->getValueByPath(hGate0, "00")
                        .approximatelyEquals(dd::COMPLEX_SQRT3_3));
    ASSERT_TRUE(dd->getValueByPath(hGate0, "10")
                        .approximatelyEquals(dd::COMPLEX_SQRT3_3));
    ASSERT_TRUE(dd->getValueByPath(hGate0, "20")
                        .approximatelyEquals(dd::COMPLEX_SQRT3_3));
    ASSERT_TRUE(dd->getValueByPath(hGate0, "30")
                        .approximatelyEquals(dd::COMPLEX_SQRT3_3));

    ASSERT_TRUE(dd->getValueByPath(hGate0, "40").approximatelyEquals(h3plus));

    ASSERT_TRUE(dd->getValueByPath(hGate0, "50").approximatelyEquals(h3minus));

    ASSERT_TRUE(dd->getValueByPath(hGate0, "60")
                        .approximatelyEquals(dd::COMPLEX_SQRT3_3));

    ASSERT_TRUE(dd->getValueByPath(hGate0, "70").approximatelyEquals(h3minus));

    ASSERT_TRUE(dd->getValueByPath(hGate0, "80").approximatelyEquals(h3plus));

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

    ASSERT_TRUE(dd->getValueByPath(hGate0, "43").approximatelyEquals(h3plus));
    ASSERT_TRUE(dd->getValueByPath(hGate0, "53").approximatelyEquals(h3minus));

    ASSERT_TRUE(dd->getValueByPath(hGate0, "63").approximatelyEquals(dd::COMPLEX_SQRT3_3));

    ASSERT_TRUE(dd->getValueByPath(hGate0, "73").approximatelyEquals(h3minus));
    ASSERT_TRUE(dd->getValueByPath(hGate0, "83").approximatelyEquals(h3plus));

    // case with higher target
    auto hGate1 = dd->makeGateDD<dd::GateMatrix>(dd::Hmat, 2, 1);
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
    auto         controlledH3 = dd->makeGateDD<dd::TritMatrix>(dd::H3(), 2, controls1, 0);

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
                        .approximatelyEquals(h3plus));
    ASSERT_TRUE(dd->getValueByPath(controlledH3, "53")
                        .approximatelyEquals(h3minus));
    ASSERT_TRUE(dd->getValueByPath(controlledH3, "63")
                        .approximatelyEquals(dd::COMPLEX_SQRT3_3));
    ASSERT_TRUE(dd->getValueByPath(controlledH3, "73")
                        .approximatelyEquals(h3minus));
    ASSERT_TRUE(dd->getValueByPath(controlledH3, "83")
                        .approximatelyEquals(h3plus));

    dd::Controls controls0{{0, 1}};
    auto         controlled0H2 =
            dd->makeGateDD<dd::GateMatrix>(dd::Hmat, 2, controls0, 1);

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

    auto id2 = dd->makeIdent(0, 1); // should be found in idTable
    EXPECT_EQ(dd->makeIdent(2), id2);

    auto id4 = dd->makeIdent(0, 3); // should use id3 and extend it
    EXPECT_EQ(dd->makeIdent(0, 3), id4);
    EXPECT_NE(table[3].nextNode, nullptr);

    auto idCached = dd->makeIdent(4);
    EXPECT_EQ(id4, idCached);
}

TEST(DDPackageTest, Multiplication) {
    auto dd =
            std::make_unique<dd::MDDPackage>(3, std::vector<std::size_t>{2, 2, 3});
    EXPECT_EQ(dd->qregisters(), 3);

    auto xGate     = dd->makeGateDD<dd::GateMatrix>(dd::Xmat, 3, 0);
    auto x3dagGate = dd->makeGateDD<dd::TritMatrix>(dd::X3dag, 3, 2);
    auto x3Gate    = dd->makeGateDD<dd::TritMatrix>(dd::X3, 3, 2);

    dd::Controls control{{0, 1}, {2, 1}};
    auto         ctrlxGate = dd->makeGateDD<dd::GateMatrix>(dd::Xmat, 3, control, 1);

    auto zeroState = dd->makeZeroState(3);

    auto evolution = dd->multiply(xGate, zeroState);

    std::cout << "\n"
              << std::endl;
    dd->printVector(evolution);

    evolution = dd->multiply(x3Gate, evolution);

    std::cout << "\n"
              << std::endl;
    dd->printVector(evolution);

    evolution = dd->multiply(ctrlxGate, evolution);

    std::cout << "\n"
              << std::endl;
    dd->printVector(evolution);

    evolution = dd->multiply(x3dagGate, evolution);

    std::cout << "\n"
              << std::endl;
    dd->printVector(evolution);

    auto basis110State = dd->makeBasisState(3, {1, 1, 0});

    ASSERT_EQ(dd->fidelity(zeroState, evolution), 0.0);
    ASSERT_EQ(dd->fidelity(evolution, basis110State), 1.0);
}

TEST(DDPackageTest, QutritBellState) {
    auto dd = std::make_unique<dd::MDDPackage>(2, std::vector<std::size_t>{3, 3});
    EXPECT_EQ(dd->qregisters(), 2);
    // Gates
    auto h3Gate = dd->makeGateDD<dd::TritMatrix>(dd::H3(), 2, 0);

    dd::Controls control01{{0, 1}};
    auto         ctrlx1Gate = dd->makeGateDD<dd::TritMatrix>(dd::X3, 2, control01, 1);

    dd::Controls control02{{0, 2}};
    auto         ctrlx2Gate = dd->makeGateDD<dd::TritMatrix>(dd::X3, 2, control02, 1);

    // Final Unitary
    auto op = dd->multiply(ctrlx1Gate, h3Gate);
    op      = dd->multiply(ctrlx2Gate, op);
    op      = dd->multiply(ctrlx2Gate, op);

    // Evolution
    auto evolution = dd->makeZeroState(2);

    std::cout << "\n"
              << std::endl;
    dd->printVector(evolution);

    evolution = dd->multiply(op, evolution);

    std::cout << "\n"
              << std::endl;
    dd->printVector(evolution);

    auto basis00State = dd->makeBasisState(2, {0, 0});
    auto basis11State = dd->makeBasisState(2, {1, 1});
    auto basis22State = dd->makeBasisState(2, {2, 2});

    ASSERT_NEAR(dd->fidelity(basis00State, evolution), 0.3333333333333333,
                dd::ComplexTable<>::tolerance());
    ASSERT_NEAR(dd->fidelity(basis11State, evolution), 0.3333333333333333,
                dd::ComplexTable<>::tolerance());
    ASSERT_NEAR(dd->fidelity(basis22State, evolution), 0.3333333333333333,
                dd::ComplexTable<>::tolerance());

    // Evolution gate times state

    evolution = dd->makeZeroState(2);

    evolution = dd->multiply(h3Gate, evolution);
    evolution = dd->multiply(ctrlx1Gate, evolution);
    evolution = dd->multiply(ctrlx2Gate, evolution);
    evolution = dd->multiply(ctrlx2Gate, evolution);

    ASSERT_NEAR(dd->fidelity(basis00State, evolution), 0.3333333333333333,
                dd::ComplexTable<>::tolerance());
    ASSERT_NEAR(dd->fidelity(basis11State, evolution), 0.3333333333333333,
                dd::ComplexTable<>::tolerance());
    ASSERT_NEAR(dd->fidelity(basis22State, evolution), 0.3333333333333333,
                dd::ComplexTable<>::tolerance());
}

TEST(DDPackageTest, GHZQutritState) {
    auto dd =
            std::make_unique<dd::MDDPackage>(3, std::vector<std::size_t>{3, 3, 3});
    EXPECT_EQ(dd->qregisters(), 3);
    // Gates
    auto h3Gate = dd->makeGateDD<dd::TritMatrix>(dd::H3(), 3, 0);

    dd::Controls control01{{0, 1}};
    auto         cX011 = dd->makeGateDD<dd::TritMatrix>(dd::X3, 3, control01, 1);
    dd::Controls control02{{0, 2}};
    auto         cX021 = dd->makeGateDD<dd::TritMatrix>(dd::X3dag, 3, control02, 1);

    dd::Controls control011{{0, 1}, {1, 1}};
    auto         cXc01l1t2 = dd->makeGateDD<dd::TritMatrix>(dd::X3, 3, control011, 2);
    dd::Controls control012{{0, 2}, {1, 2}};
    auto         cX0122 = dd->makeGateDD<dd::TritMatrix>(dd::X3dag, 3, control012, 2);

    //auto testvec = dd->getVectorizedMatrix(cX0122);

    // Evolution
    auto evolution = dd->makeZeroState(3);

    std::cout << "\n"
              << std::endl;
    dd->printVector(evolution);

    evolution = dd->multiply(h3Gate, evolution);

    evolution = dd->multiply(cX011, evolution);
    evolution = dd->multiply(cX021, evolution);

    std::cout << "\n"
              << std::endl;
    dd->printVector(evolution);

    evolution = dd->multiply(cXc01l1t2, evolution);
    std::cout << "\n"
              << std::endl;
    dd->printVector(evolution);

    evolution = dd->multiply(cX0122, evolution);

    std::cout << "\n"
              << std::endl;
    dd->printVector(evolution);

    auto basis00State = dd->makeBasisState(3, {0, 0, 0});
    auto basis11State = dd->makeBasisState(3, {1, 1, 1});
    auto basis22State = dd->makeBasisState(3, {2, 2, 2});

    ASSERT_NEAR(dd->fidelity(basis00State, evolution), 0.3333333333333333,
                dd::ComplexTable<>::tolerance());
    ASSERT_NEAR(dd->fidelity(basis11State, evolution), 0.3333333333333333,
                dd::ComplexTable<>::tolerance());
    ASSERT_NEAR(dd->fidelity(basis22State, evolution), 0.3333333333333333,
                dd::ComplexTable<>::tolerance());
}

TEST(DDPackageTest, GHZQutritStateScaled) {
    for (auto i = 20U; i < 21U; i++) {
        auto                     test = 0;
        std::vector<std::size_t> init(i, 3);
        auto                     dd = std::make_unique<dd::MDDPackage>(i, init);
        EXPECT_EQ(dd->qregisters(), i);

        // Gates
        auto                               h3Gate = dd->makeGateDD<dd::TritMatrix>(dd::H3(), i, 0);
        std::vector<dd::MDDPackage::mEdge> gates  = {};

        for (auto target = 1U; target < i; target++) {
            dd::Controls target1{};
            dd::Controls target2{};

            for (auto control = 0U; control < target; control++) {
                const dd::Control c1{static_cast<dd::QuantumRegister>(control), 1};
                const dd::Control c2{static_cast<dd::QuantumRegister>(control), 2};
                target1.insert(c1);
                target2.insert(c2);
            }

            gates.push_back(
                    dd->makeGateDD<dd::TritMatrix>(dd::X3, i, target1, target));
            gates.push_back(
                    dd->makeGateDD<dd::TritMatrix>(dd::X3dag, i, target2, target));
        }

        auto evolution = dd->makeZeroState(i);
        evolution      = dd->multiply(h3Gate, evolution);

        for (auto& gate: gates) {
            evolution = dd->multiply(gate, evolution);
        }
        if (dd->qregisters() < 10) {
            std::cout << "\n"
                      << std::endl;
            dd->printVector(evolution);
        }

        auto basis00State = dd->makeBasisState(i, std::vector<size_t>(i, 0));
        auto basis11State = dd->makeBasisState(i, std::vector<size_t>(i, 1));
        auto basis22State = dd->makeBasisState(i, std::vector<size_t>(i, 2));

        ASSERT_NEAR(dd->fidelity(basis00State, evolution), 0.3333333333333333,
                    dd::ComplexTable<>::tolerance());
        ASSERT_NEAR(dd->fidelity(basis11State, evolution), 0.3333333333333333,
                    dd::ComplexTable<>::tolerance());
        ASSERT_NEAR(dd->fidelity(basis22State, evolution), 0.3333333333333333,
                    dd::ComplexTable<>::tolerance());
    }
}

TEST(DDPackageTest, RandomCircuits) {
    int          width = 6;
    unsigned int depth = 1L;
    size_t       maxD  = 5;

    std::random_device rd;        // obtain a random number from hardware
    std::mt19937       gen(rd()); // seed the generator

    std::vector<std::size_t> particles = {};

    std::uniform_int_distribution<> dimdistr(2, maxD);
    for (auto i = 0; i < width; i++) {
        particles.push_back(dimdistr(gen));
    }

    auto dd = std::make_unique<dd::MDDPackage>(width, particles);

    EXPECT_EQ(dd->qregisters(), width);

    auto evolution = dd->makeZeroState(width);

    std::uniform_int_distribution<>  pickbool(0, 1);
    std::uniform_int_distribution<>  pickcontrols(1, width - 1);
    std::uniform_real_distribution<> angles(0.0, 2. * dd::PI);

    for (auto timeStep = 0U; timeStep < depth; timeStep++) {
        for (int line = 0; line < width; line++) {
            //chose if local gate or entangling gate
            auto randomChoice = pickbool(gen);

            if (randomChoice == 0) { //local op

                auto localChoice = pickbool(gen);

                if (localChoice == 0) { //hadamard
                    if (particles.at(line) == 2) {
                        auto chosenGate = dd->makeGateDD<dd::GateMatrix>(dd::H(), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 3) {
                        auto chosenGate = dd->makeGateDD<dd::TritMatrix>(dd::H3(), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 4) {
                        auto chosenGate = dd->makeGateDD<dd::QuartMatrix>(dd::H4(), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 5) {
                        auto chosenGate = dd->makeGateDD<dd::QuintMatrix>(dd::H5(), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    }
                } else { //givens

                    if (particles.at(line) == 2) {
                        long double theta      = 0.;
                        long double phi        = 0.;
                        auto        chosenGate = dd->makeGateDD<dd::GateMatrix>(dd::RXY(theta, phi), width, line);
                        evolution              = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 3) {
                        long double                     theta = angles(gen);
                        long double                     phi   = angles(gen);
                        std::uniform_int_distribution<> picklevel(0, 2);

                        int levelA = picklevel(gen);
                        int levelB = (levelA + 1) % 3;
                        if (levelA > levelB) {
                            auto temp = levelA;
                            levelA    = levelB;
                            levelB    = temp;
                        }

                        auto chosenGate = dd->makeGateDD<dd::TritMatrix>(dd::RXY3(theta, phi, levelA, levelB), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 4) {
                        long double                     theta = angles(gen);
                        long double                     phi   = angles(gen);
                        std::uniform_int_distribution<> picklevel(0, 3);

                        int levelA = picklevel(gen);
                        int levelB = (levelA + 1) % 4;
                        if (levelA > levelB) {
                            auto temp = levelA;
                            levelA    = levelB;
                            levelB    = temp;
                        }

                        auto chosenGate = dd->makeGateDD<dd::QuartMatrix>(dd::RXY4(theta, phi, levelA, levelB), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 5) {
                        long double                     theta = angles(gen);
                        long double                     phi   = angles(gen);
                        std::uniform_int_distribution<> picklevel(0, 4);

                        int levelA = picklevel(gen);
                        int levelB = (levelA + 1) % 5;
                        if (levelA > levelB) {
                            auto temp = levelA;
                            levelA    = levelB;
                            levelB    = temp;
                        }

                        auto chosenGate = dd->makeGateDD<dd::QuintMatrix>(dd::RXY5(theta, phi, levelA, levelB), width, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    }
                }
            } else { //entangling

                auto entChoice        = pickbool(gen);
                auto numberOfControls = pickcontrols(gen);

            std:
                std::vector<int> controlLines;

                for (auto i = 0U; i < width; i++) {
                    if (i != line) controlLines.push_back(i);
                }

                std::shuffle(begin(controlLines), end(controlLines), gen);
                std::vector<int> controlParticles(controlLines.begin(), controlLines.begin() + numberOfControls);
                std::sort(controlParticles.begin(), controlParticles.end());

                dd::Controls control{};
                for (auto i = 0U; i < numberOfControls; i++) {
                    std::uniform_int_distribution<> picklevel(0, particles.at(controlParticles.at(i)) - 1);
                    auto                            level = picklevel(gen);

                    const dd::Control c{static_cast<dd::QuantumRegister>(controlParticles.at(i)), static_cast<dd::Control::Type>(level)};
                    control.insert(c);
                }

                if (entChoice == 0) { // CEX based
                    //selection of controls

                    if (particles.at(line) == 2) {
                        long double theta      = angles(gen);
                        long double phi        = angles(gen);
                        auto        chosenGate = dd->makeGateDD<dd::GateMatrix>(dd::RXY(theta, phi), width, control, line);
                        evolution              = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 3) {
                        std::uniform_int_distribution<> picklevel(0, 2);

                        long double theta  = angles(gen);
                        long double phi    = angles(gen);
                        int         levelA = picklevel(gen);
                        int         levelB = (levelA + 1) % 3;
                        if (levelA > levelB) {
                            auto temp = levelA;
                            levelA    = levelB;
                            levelB    = temp;
                        }

                        auto chosenGate = dd->makeGateDD<dd::TritMatrix>(dd::RXY3(theta, phi, levelA, levelB), width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 4) {
                        std::uniform_int_distribution<> picklevel(0, 3);

                        long double theta  = angles(gen);
                        long double phi    = angles(gen);
                        int         levelA = picklevel(gen);
                        int         levelB = (levelA + 1) % 4;
                        if (levelA > levelB) {
                            auto temp = levelA;
                            levelA    = levelB;
                            levelB    = temp;
                        }

                        auto chosenGate = dd->makeGateDD<dd::QuartMatrix>(dd::RXY4(theta, phi, levelA, levelB), width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 5) {
                        std::uniform_int_distribution<> picklevel(0, 4);

                        long double theta  = angles(gen);
                        long double phi    = angles(gen);
                        int         levelA = picklevel(gen);
                        int         levelB = (levelA + 1) % 5;
                        if (levelA > levelB) {
                            auto temp = levelA;
                            levelA    = levelB;
                            levelB    = temp;
                        }

                        auto chosenGate = dd->makeGateDD<dd::QuintMatrix>(dd::RXY5(theta, phi, levelA, levelB), width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    }
                } else { // Controlled clifford
                    if (particles.at(line) == 2) {
                        auto chosenGate = dd->makeGateDD<dd::GateMatrix>(dd::Xmat, width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 3) {
                        auto chosenGate = dd->makeGateDD<dd::TritMatrix>(dd::X3, width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 4) {
                        auto chosenGate = dd->makeGateDD<dd::QuartMatrix>(dd::X4, width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    } else if (particles.at(line) == 5) {
                        auto chosenGate = dd->makeGateDD<dd::QuintMatrix>(dd::X5, width, control, line);
                        evolution       = dd->multiply(chosenGate, evolution);
                    }
                }
            }
        }
    }
    auto basis0State = dd->makeBasisState(width, std::vector<size_t>(width, 0));
    dd->printVector(evolution);
    //ASSERT_NEAR(dd->fidelity(basis0State, evolution), 0.3333333333333333, dd::ComplexTable<>::tolerance());
}
