/*
 * This file is part of the MQT DD Package which is released under the MIT
 * license. See file README.md or go to
 * http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DD_PACKAGE_GATEMATRIXDEFINITIONS_H
#define DD_PACKAGE_GATEMATRIXDEFINITIONS_H

#include "ComplexValue.hpp"
#include "Definitions.hpp"

#include <cmath>
#include <vector>

namespace dd {
    // Complex constants
    constexpr ComplexValue COMPLEX_ONE       = {1., 0.};
    constexpr ComplexValue COMPLEX_MONE      = {-1., 0.};
    constexpr ComplexValue COMPLEX_ZERO      = {0., 0.};
    constexpr ComplexValue COMPLEX_I         = {0., 1.};
    constexpr ComplexValue COMPLEX_MI        = {0., -1.};
    constexpr ComplexValue COMPLEX_SQRT2_2   = {SQRT2_2, 0.};
    constexpr ComplexValue COMPLEX_MSQRT2_2  = {-SQRT2_2, 0.};
    constexpr ComplexValue COMPLEX_ISQRT2_2  = {0., SQRT2_2};
    constexpr ComplexValue COMPLEX_MISQRT2_2 = {0., -SQRT2_2};
    constexpr ComplexValue COMPLEX_1PLUSI    = {SQRT2_2, SQRT2_2};
    constexpr ComplexValue COMPLEX_1MINUSI   = {SQRT2_2, -SQRT2_2};
    constexpr ComplexValue COMPLEX_1PLUSI_2  = {0.5, 0.5};
    constexpr ComplexValue COMPLEX_1MINUSI_2 = {0.5, -0.5};

    using GateMatrix = std::array<ComplexValue, EDGE2>;
    using TritMatrix = std::array<ComplexValue, EDGE3>;

    constexpr GateMatrix Imat{COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE};
    constexpr GateMatrix Hmat{COMPLEX_SQRT2_2, COMPLEX_SQRT2_2, COMPLEX_SQRT2_2,
                              COMPLEX_MSQRT2_2};
    constexpr GateMatrix Xmat{COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ONE, COMPLEX_ZERO};
    constexpr GateMatrix Ymat{COMPLEX_ZERO, COMPLEX_MI, COMPLEX_I, COMPLEX_ZERO};
    constexpr GateMatrix Zmat{COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
                              COMPLEX_MONE};
    constexpr GateMatrix Smat{COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_I};
    constexpr GateMatrix Sdagmat{COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
                                 COMPLEX_MI};
    constexpr GateMatrix Tmat{COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
                              COMPLEX_1PLUSI};
    constexpr GateMatrix Tdagmat{COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
                                 COMPLEX_1MINUSI};
    constexpr GateMatrix SXmat{COMPLEX_1PLUSI_2, COMPLEX_1MINUSI_2,
                               COMPLEX_1MINUSI_2, COMPLEX_1PLUSI_2};
    constexpr GateMatrix SXdagmat{COMPLEX_1MINUSI_2, COMPLEX_1PLUSI_2,
                                  COMPLEX_1PLUSI_2, COMPLEX_1MINUSI_2};
    constexpr GateMatrix Vmat{COMPLEX_SQRT2_2, COMPLEX_MISQRT2_2, COMPLEX_MISQRT2_2,
                              COMPLEX_SQRT2_2};
    constexpr GateMatrix Vdagmat{COMPLEX_SQRT2_2, COMPLEX_ISQRT2_2,
                                 COMPLEX_ISQRT2_2, COMPLEX_SQRT2_2};

    // TODO FIX THE CONSTANTS
    constexpr ComplexValue COMPLEX_SQRT3_3 = {0.57735026918962576450914878050195745565, 0.};
    constexpr ComplexValue COMPLEX_SQRT3_3 = {0.57735026918962576450914878050195745565, 0.};
    constexpr ComplexValue COMPLEX_SQRT3_3 = {0.57735026918962576450914878050195745565, 0.};
    constexpr ComplexValue COMPLEX_SQRT3_3 = {0.57735026918962576450914878050195745565, 0.};

    constexpr TritMatrix H3{COMPLEX_SQRT3_3,
                            COMPLEX_SQRT3_3,
                            COMPLEX_SQRT3_3,
                            COMPLEX_SQRT3_3,
                            dd::ComplexValue{-0.28867513, 0.5},
                            dd::ComplexValue{-0.28867513, -0.5},
                            COMPLEX_SQRT3_3,
                            dd::ComplexValue{-0.28867513, -0.5},
                            dd::ComplexValue{-0.28867513, 0.5}};

    constexpr TritMatrix X3dag{COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO,
                               COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE,
                               COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO};

    constexpr TritMatrix X3{COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE,
                            COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
                            COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO};

    constexpr TritMatrix PI02{COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_MONE,
                              COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO,
                              COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO};
    constexpr TritMatrix P_I02DAG{COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE,
                                  COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO,
                                  COMPLEX_MONE, COMPLEX_ZERO, COMPLEX_ZERO};
    constexpr TritMatrix X01{COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO,
                             COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
                             COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE};
    constexpr TritMatrix Z01{COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
                             COMPLEX_ZERO, COMPLEX_MONE, COMPLEX_ZERO,
                             COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE};

    inline TritMatrix U301(fp lambda, fp phi, fp theta) {
        return TritMatrix{{{std::cos(theta / 2.), 0.},
                           {-std::cos(lambda) * std::sin(theta / 2.),
                            -std::sin(lambda) * std::sin(theta / 2.)},
                           COMPLEX_ZERO,
                           {std::cos(phi) * std::sin(theta / 2.),
                            std::sin(phi) * std::sin(theta / 2.)},
                           {std::cos(lambda + phi) * std::cos(theta / 2.),
                            std::sin(lambda + phi) * std::cos(theta / 2.)},
                           COMPLEX_ZERO,
                           COMPLEX_ZERO,
                           COMPLEX_ZERO,
                           COMPLEX_ONE}};
    }

    // TODO CHECK DEFINITIONS ARE WEIRD
    inline GateMatrix U3mat(fp lambda, fp phi, fp theta) {
        return GateMatrix{{{std::cos(theta / 2.), 0.},
                           {-std::cos(lambda) * std::sin(theta / 2.),
                            -std::sin(lambda) * std::sin(theta / 2.)},
                           {std::cos(phi) * std::sin(theta / 2.),
                            std::sin(phi) * std::sin(theta / 2.)},
                           {std::cos(lambda + phi) * std::cos(theta / 2.),
                            std::sin(lambda + phi) * std::cos(theta / 2.)}}};
    }

    inline GateMatrix U2mat(fp lambda, fp phi) {
        return GateMatrix{
                COMPLEX_SQRT2_2,
                {-std::cos(lambda) * SQRT2_2, -std::sin(lambda) * SQRT2_2},
                {std::cos(phi) * SQRT2_2, std::sin(phi) * SQRT2_2},
                {std::cos(lambda + phi) * SQRT2_2, std::sin(lambda + phi) * SQRT2_2}};
    }

    inline GateMatrix Phasemat(fp lambda) {
        return GateMatrix{COMPLEX_ONE,
                          COMPLEX_ZERO,
                          COMPLEX_ZERO,
                          {std::cos(lambda), std::sin(lambda)}};
    }

    inline GateMatrix RXmat(fp lambda) {
        return GateMatrix{{{std::cos(lambda / 2.), 0.},
                           {0., -std::sin(lambda / 2.)},
                           {0., -std::sin(lambda / 2.)},
                           {std::cos(lambda / 2.), 0.}}};
    }

    inline GateMatrix RYmat(fp lambda) {
        return GateMatrix{{{std::cos(lambda / 2.), 0.},
                           {-std::sin(lambda / 2.), 0.},
                           {std::sin(lambda / 2.), 0.},
                           {std::cos(lambda / 2.), 0.}}};
    }

    inline GateMatrix RZmat(fp lambda) {
        return GateMatrix{{{std::cos(lambda / 2.), -std::sin(lambda / 2.)},
                           COMPLEX_ZERO,
                           COMPLEX_ZERO,
                           {std::cos(lambda / 2.), std::sin(lambda / 2.)}}};
    }
} // namespace dd
#endif // DD_PACKAGE_GATEMATRIXDEFINITIONS_H
