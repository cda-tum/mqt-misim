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
#include <stdexcept>
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

    constexpr ComplexValue COMPLEX_SQRT3_3 = {SQRT3_3, 0.};
    constexpr ComplexValue COMPLEX_SQRT4_4 = {SQRT4_4, 0.};
    constexpr ComplexValue COMPLEX_SQRT5_5 = {SQRT5_5, 0.};

    using GateMatrix  = std::array<ComplexValue, EDGE2>;
    using TritMatrix  = std::array<ComplexValue, EDGE3>;
    using QuartMatrix = std::array<ComplexValue, EDGE4>;
    using QuintMatrix = std::array<ComplexValue, EDGE5>;

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

    inline GateMatrix Pimat(size_t i) {
        GateMatrix zero    = {COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO};
        zero.at(i + i * 2) = COMPLEX_ONE;
        return zero;
    }

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

    inline GateMatrix H() {
        return GateMatrix{COMPLEX_SQRT2_2,
                          COMPLEX_SQRT2_2,
                          COMPLEX_SQRT2_2,
                          COMPLEX_MSQRT2_2};
    }
    inline GateMatrix RXY(fp theta, fp phi) {
        GateMatrix rotation = {
                dd::ComplexValue{std::cos(theta / 2.), 0.},
                dd::ComplexValue{-std::sin(theta / 2.) * std::sin(phi), -std::sin(theta / 2.) * std::cos(phi)},
                dd::ComplexValue{std::sin(theta / 2.) * std::sin(phi), -std::sin(theta / 2.) * std::cos(phi)},
                dd::ComplexValue{std::cos(theta / 2.), 0.}};
        return rotation;
    }

    constexpr TritMatrix I3{COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
                            COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO,
                            COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE};

    inline TritMatrix Pi3(size_t i) {
        TritMatrix zero    = {COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
                              COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
                              COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO};
        zero.at(i + i * 3) = COMPLEX_ONE;
        return zero;
    }

    inline TritMatrix H3() {
        return TritMatrix{COMPLEX_SQRT3_3,
                          COMPLEX_SQRT3_3,
                          COMPLEX_SQRT3_3,

                          COMPLEX_SQRT3_3,
                          COMPLEX_SQRT3_3 * dd::ComplexValue{std::cos(2. * PI / 3.), std::sin(2. * PI / 3.)},
                          COMPLEX_SQRT3_3 * dd::ComplexValue{std::cos(4. * PI / 3.), std::sin(4. * PI / 3.)},

                          COMPLEX_SQRT3_3,
                          COMPLEX_SQRT3_3 * dd::ComplexValue{std::cos(4. * PI / 3.), std::sin(4. * PI / 3.)},
                          COMPLEX_SQRT3_3 * dd::ComplexValue{std::cos(2. * PI / 3.), std::sin(2. * PI / 3.)}};
    }

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

    inline TritMatrix RXY3(fp theta, fp phi, size_t leva, size_t levb) {
        if (leva > levb or leva > 2 or levb > 2) {
            throw std::invalid_argument("LEV A cannot be higher than  LEV B");
        }
        TritMatrix identity = {
                COMPLEX_ONE,
                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ONE,
                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ONE};
        identity.at(3 * leva + leva) = dd::ComplexValue{std::cos(theta / 2.), 0.};
        identity.at(3 * leva + levb) = dd::ComplexValue{-std::sin(theta / 2.) * std::sin(phi), -std::sin(theta / 2.) * std::cos(phi)};
        identity.at(3 * levb + leva) = dd::ComplexValue{std::sin(theta / 2.) * std::sin(phi), -std::sin(theta / 2.) * std::cos(phi)};
        identity.at(3 * levb + levb) = dd::ComplexValue{std::cos(theta / 2.), 0.};
        return identity;
    }

    constexpr QuartMatrix I4{
            COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
            COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
            COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO,
            COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE};

    inline QuartMatrix Pi4(size_t i) {
        QuartMatrix zero   = {COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
                              COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
                              COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
                              COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO};
        zero.at(i + i * 4) = COMPLEX_ONE;
        return zero;
    }

    inline QuartMatrix H4() {
        return QuartMatrix{COMPLEX_SQRT4_4,
                           COMPLEX_SQRT4_4,
                           COMPLEX_SQRT4_4,
                           COMPLEX_SQRT4_4,

                           COMPLEX_SQRT4_4,
                           COMPLEX_SQRT4_4 * dd::ComplexValue{std::cos(2. * PI / 4.), std::sin(2. * PI / 4.)},
                           COMPLEX_SQRT4_4 * dd::ComplexValue{std::cos(2. * 2. * PI / 4.), std::sin(2. * 2. * PI / 4.)},
                           COMPLEX_SQRT4_4 * dd::ComplexValue{std::cos(3. * 2. * PI / 4.), std::sin(3. * 2. * PI / 4.)},

                           COMPLEX_SQRT4_4,
                           COMPLEX_SQRT4_4 * dd::ComplexValue{std::cos(2. * 2. * PI / 4.), std::sin(2. * 2. * PI / 4.)},
                           COMPLEX_SQRT4_4 * dd::ComplexValue{std::cos(4. * 2. * PI / 4.), std::sin(4. * 2. * PI / 4.)},
                           COMPLEX_SQRT4_4 * dd::ComplexValue{std::cos(6. * 2. * PI / 4.), std::sin(6. * 2. * PI / 4.)},

                           COMPLEX_SQRT4_4,
                           COMPLEX_SQRT4_4 * dd::ComplexValue{std::cos(3. * 2. * PI / 4.), std::sin(3. * 2. * PI / 4.)},
                           COMPLEX_SQRT4_4 * dd::ComplexValue{std::cos(6. * 2. * PI / 4.), std::sin(6. * 2. * PI / 4.)},
                           COMPLEX_SQRT4_4 * dd::ComplexValue{std::cos(9. * 2. * PI / 4.), std::sin(9. * 2. * PI / 4.)}};
    }
    inline QuartMatrix RXY4(fp theta, fp phi, size_t leva, size_t levb) {
        if (leva > levb or leva > 3 or levb > 3) {
            throw std::invalid_argument("LEV A cannot be higher than  LEV B");
        }
        QuartMatrix identity = {
                COMPLEX_ONE,
                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ZERO,

                COMPLEX_ZERO,
                COMPLEX_ONE,
                COMPLEX_ZERO,
                COMPLEX_ZERO,

                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ONE,
                COMPLEX_ZERO,

                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ONE};
        identity.at(4 * leva + leva) = dd::ComplexValue{std::cos(theta / 2.), 0.};
        identity.at(4 * leva + levb) = dd::ComplexValue{-std::sin(theta / 2.) * std::sin(phi), -std::sin(theta / 2.) * std::cos(phi)};
        identity.at(4 * levb + leva) = dd::ComplexValue{std::sin(theta / 2.) * std::sin(phi), -std::sin(theta / 2.) * std::cos(phi)};
        identity.at(4 * levb + levb) = dd::ComplexValue{std::cos(theta / 2.), 0.};
        return identity;
    }
    constexpr QuartMatrix X4dag{COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
                                COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO,
                                COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE,
                                COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO};

    constexpr QuartMatrix X4{COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE,
                             COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
                             COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
                             COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO};

    constexpr QuintMatrix I5{
            COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
            COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
            COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
            COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO,
            COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE};

    inline QuintMatrix Pi5(size_t i) {
        QuintMatrix zero   = {COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
                              COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
                              COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
                              COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
                              COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO};
        zero.at(i + i * 5) = COMPLEX_ONE;
        return zero;
    }
    inline QuintMatrix H5() {
        return QuintMatrix{COMPLEX_SQRT5_5,
                           COMPLEX_SQRT5_5,
                           COMPLEX_SQRT5_5,
                           COMPLEX_SQRT5_5,
                           COMPLEX_SQRT5_5,

                           COMPLEX_SQRT5_5,
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(2. * PI / 5.), std::sin(2. * PI / 5.)},
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(2. * 2. * PI / 5.), std::sin(2. * 2. * PI / 5.)},
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(3. * 2. * PI / 5.), std::sin(3. * 2. * PI / 5.)},
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(4. * 2. * PI / 5.), std::sin(4. * 2. * PI / 5.)},

                           COMPLEX_SQRT5_5,
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(2. * 2. * PI / 5.), std::sin(2. * 2. * PI / 5.)},
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(4. * 2. * PI / 5.), std::sin(4. * 2. * PI / 5.)},
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(6. * 2. * PI / 5.), std::sin(6. * 2. * PI / 5.)},
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(8. * 2. * PI / 5.), std::sin(8. * 2. * PI / 5.)},

                           COMPLEX_SQRT5_5,
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(3. * 2. * PI / 5.), std::sin(3. * 2. * PI / 5.)},
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(6. * 2. * PI / 5.), std::sin(6. * 2. * PI / 5.)},
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(9. * 2. * PI / 5.), std::sin(9. * 2. * PI / 5.)},
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(12. * 2. * PI / 5.), std::sin(12. * 2. * PI / 5.)},

                           COMPLEX_SQRT5_5,
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(4. * 2. * PI / 5.), std::sin(4. * 2. * PI / 5.)},
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(8. * 2. * PI / 5.), std::sin(8. * 2. * PI / 5.)},
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(12. * 2. * PI / 5.), std::sin(12. * 2. * PI / 5.)},
                           COMPLEX_SQRT5_5 * dd::ComplexValue{std::cos(16. * 2. * PI / 5.), std::sin(16. * 2. * PI / 5.)}};
    }

    inline QuintMatrix RXY5(fp theta, fp phi, size_t leva, size_t levb) {
        if (leva > levb or leva > 4 or levb > 4) {
            throw std::invalid_argument("LEV A cannot be higher than  LEV B");
        }
        QuintMatrix identity = {
                COMPLEX_ONE,
                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ZERO,

                COMPLEX_ZERO,
                COMPLEX_ONE,
                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ZERO,

                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ONE,
                COMPLEX_ZERO,
                COMPLEX_ZERO,

                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ONE,
                COMPLEX_ZERO,

                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ZERO,
                COMPLEX_ONE};
        identity.at(5 * leva + leva) = dd::ComplexValue{std::cos(theta / 2.), 0.};
        identity.at(5 * leva + levb) = dd::ComplexValue{-std::sin(theta / 2.) * std::sin(phi), -std::sin(theta / 2.) * std::cos(phi)};
        identity.at(5 * levb + leva) = dd::ComplexValue{std::sin(theta / 2.) * std::sin(phi), -std::sin(theta / 2.) * std::cos(phi)};
        identity.at(5 * levb + levb) = dd::ComplexValue{std::cos(theta / 2.), 0.};
        return identity;
    }

    constexpr QuintMatrix X5dag{
            COMPLEX_ZERO,
            COMPLEX_ONE,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ONE,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ONE,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ONE,
            COMPLEX_ONE,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
            COMPLEX_ZERO,
    };

    constexpr QuintMatrix X5{COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE,
                             COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
                             COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO,
                             COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO, COMPLEX_ZERO,
                             COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ZERO, COMPLEX_ONE, COMPLEX_ZERO};
} // namespace dd
#endif // DD_PACKAGE_GATEMATRIXDEFINITIONS_H
