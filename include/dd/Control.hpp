/*
 * This file is part of the MQT DD Package which is released under the MIT license.
 * See file README.md or go to http://iic.jku.at/eda/research/quantum_dd/ for more information.
 */

#ifndef DD_PACKAGE_CONTROL_HPP
#define DD_PACKAGE_CONTROL_HPP

#include "Definitions.hpp"

#include <set>

namespace dd {
    struct Control {
        // make dd -> op you need -> polarized control, unsigned integer
        using Type = std::uint_fast8_t;

        QuantumRegister quantum_register{};
        Type  type = 1; // default value of a control level
    };

    inline bool operator<(const Control& lhs, const Control& rhs) {
        return lhs.quantum_register < rhs.quantum_register || (lhs.quantum_register == rhs.quantum_register && lhs.type < rhs.type);
    }

    inline bool operator==(const Control& lhs, const Control& rhs) {
        return lhs.quantum_register == rhs.quantum_register && lhs.type == rhs.type;
    }

    inline bool operator!=(const Control& lhs, const Control& rhs) {
        return !(lhs == rhs);
    }

    // this allows a set of controls to be indexed by a quantum register, namely a qudit w/ d>=2
    struct CompareControl {
        using is_transparent = void;

        inline bool operator()(const Control& lhs, const Control& rhs) const {
            return lhs < rhs;
        }

        inline bool operator()(QuantumRegister lhs, const Control& rhs) const {
            return lhs < rhs.quantum_register;
        }

        inline bool operator()(const Control& lhs, QuantumRegister rhs) const {
            return lhs.quantum_register < rhs;
        }
    };
    using Controls = std::set<Control, CompareControl>;

    inline namespace literals {
        inline Control operator""_pc(unsigned long long int qreg) { return {static_cast<QuantumRegister>(qreg)}; }
        //inline Control operator""_nc(unsigned long long int qreg) { return {static_cast<QuantumRegister>(qreg), Control::Type::neg}; }
        //inline Control operator""_nc(unsigned long long int qreg) { return {static_cast<QuantumRegister>(qreg), Control::Type::neg}; }
    } // namespace literals
} // namespace dd

#endif //DD_PACKAGE_CONTROL_HPP
