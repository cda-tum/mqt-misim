/*
 * This file is part of the MQT DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DD_PACKAGE_GATEMATRIXDEFINITIONS_H
#define DD_PACKAGE_GATEMATRIXDEFINITIONS_H

#include "ComplexValue.hpp"
#include "Definitions.hpp"

#include <vector>
#include <cmath>

namespace dd {
    // Complex constants
    constexpr ComplexValue complex_one       = {1., 0.};
    constexpr ComplexValue complex_mone      = {-1., 0.};
    constexpr ComplexValue complex_zero      = {0., 0.};
    constexpr ComplexValue complex_i         = {0., 1.};
    constexpr ComplexValue complex_mi        = {0., -1.};
    constexpr ComplexValue complex_SQRT2_2   = {SQRT2_2, 0.};
    constexpr ComplexValue complex_mSQRT2_2  = {-SQRT2_2, 0.};
    constexpr ComplexValue complex_iSQRT2_2  = {0., SQRT2_2};
    constexpr ComplexValue complex_miSQRT2_2 = {0., -SQRT2_2};
    constexpr ComplexValue complex_1plusi    = {SQRT2_2, SQRT2_2};
    constexpr ComplexValue complex_1minusi   = {SQRT2_2, -SQRT2_2};
    constexpr ComplexValue complex_1plusi_2  = {0.5, 0.5};
    constexpr ComplexValue complex_1minusi_2 = {0.5, -0.5};


    using GateMatrix = std::array<ComplexValue, EDGE2>;
    using TritMatrix = std::array<ComplexValue, EDGE3>;

    constexpr TritMatrix X{complex_zero, complex_zero, complex_one, complex_zero, complex_one, complex_zero, complex_one, complex_zero, complex_zero};
    constexpr TritMatrix PI02{complex_zero, complex_zero, complex_mone, complex_zero, complex_one, complex_zero, complex_one, complex_zero, complex_zero};
    constexpr TritMatrix PI02dag{complex_zero, complex_zero, complex_one, complex_zero, complex_one, complex_zero, complex_mone, complex_zero, complex_zero};
    constexpr TritMatrix X01{complex_zero, complex_one, complex_zero, complex_one, complex_zero, complex_zero, complex_zero, complex_zero, complex_one};
    constexpr TritMatrix Z01{complex_one, complex_zero, complex_zero,  complex_zero, complex_mone, complex_zero, complex_zero, complex_zero, complex_one};

    inline TritMatrix U3embedded(fp lambda, fp phi, fp theta) {
        return TritMatrix{
                {
                    {std::cos(theta / 2.), 0.},
                    {-std::cos(lambda) * std::sin(theta / 2.), -std::sin(lambda) * std::sin(theta / 2.)},
                         complex_zero,
                    {std::cos(phi) * std::sin(theta / 2.), std::sin(phi) * std::sin(theta / 2.)},
                    {std::cos(lambda + phi) * std::cos(theta / 2.), std::sin(lambda + phi) * std::cos(theta / 2.)},
                         complex_zero,
                         complex_zero, complex_zero, complex_one
                         }
                         };
    }

    inline GateMatrix U3mat(fp lambda, fp phi, fp theta) {
        return GateMatrix{{{std::cos(theta / 2.), 0.},
                           {-std::cos(lambda) * std::sin(theta / 2.), -std::sin(lambda) * std::sin(theta / 2.)},
                           {std::cos(phi) * std::sin(theta / 2.), std::sin(phi) * std::sin(theta / 2.)},
                           {std::cos(lambda + phi) * std::cos(theta / 2.), std::sin(lambda + phi) * std::cos(theta / 2.)}}};
    }

    inline GateMatrix U2mat(fp lambda, fp phi) {
        return GateMatrix{complex_SQRT2_2,
                          {-std::cos(lambda) * SQRT2_2, -std::sin(lambda) * SQRT2_2},
                          {std::cos(phi) * SQRT2_2, std::sin(phi) * SQRT2_2},
                          {std::cos(lambda + phi) * SQRT2_2, std::sin(lambda + phi) * SQRT2_2}};
    }

    inline GateMatrix Phasemat(fp lambda) {
        return GateMatrix{complex_one, complex_zero, complex_zero, {std::cos(lambda), std::sin(lambda)}};
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
                           complex_zero,
                           complex_zero,
                           {std::cos(lambda / 2.), std::sin(lambda / 2.)}}};
    }
} // namespace dd
#endif //DD_PACKAGE_GATEMATRIXDEFINITIONS_H
