/*
 * This file is part of the MQT DD Package which is released under the MIT license.
 * See file README.md or go to https://www.cda.cit.tum.de/research/quantum_dd/ for more information.
 */

#ifndef DD_PACKAGE_COMPLEX_HPP
#define DD_PACKAGE_COMPLEX_HPP

#include "ComplexTable.hpp"
#include "ComplexValue.hpp"

#include <cstddef>
#include <iostream>
#include <utility>

namespace dd {
    using CTEntry = ComplexTable<>::Entry;

    struct Complex {
        CTEntry* real;
        CTEntry* img;

        static Complex zero;
        static Complex one;

        void setVal(const Complex& complex_num) const {
            real->value = CTEntry::val(complex_num.real);
            img->value  = CTEntry::val(complex_num.img);
        }

        [[nodiscard]] inline bool approximatelyEquals(const Complex& complex_num) const {
            return CTEntry::approximatelyEquals(real, complex_num.real) && CTEntry::approximatelyEquals(img, complex_num.img);
        };

        [[nodiscard]] inline bool approximatelyZero() const {
            return CTEntry::approximatelyZero(real) && CTEntry::approximatelyZero(img);
        }

        [[nodiscard]] inline bool approximatelyOne() const {
            return CTEntry::approximatelyOne(real) && CTEntry::approximatelyZero(img);
        }

        inline bool operator==(const Complex& other) const {
            return real == other.real && img == other.img;
        }

        inline bool operator!=(const Complex& other) const {
            return !operator==(other);
        }

        [[nodiscard]] std::string toString(bool formatted = true, int precision = -1) const {
            return ComplexValue::toString(CTEntry::val(real), CTEntry::val(img), formatted, precision);
        }

        void writeBinary(std::ostream& os) const {
            CTEntry::writeBinary(real, os);
            CTEntry::writeBinary(img, os);
        }
    };

    inline std::ostream& operator<<(std::ostream& os, const Complex& complex_num) {
        return os << complex_num.toString();
    }

    inline Complex Complex::zero{&ComplexTable<>::zero, &ComplexTable<>::zero};
    inline Complex Complex::one{&ComplexTable<>::one, &ComplexTable<>::zero};
} // namespace dd

namespace std {
    template<>
    struct hash<dd::Complex> {
        std::size_t operator()(dd::Complex const& complex_num) const noexcept {
            auto h1 = dd::murmur64(reinterpret_cast<std::size_t>(complex_num.real));
            auto h2 = dd::murmur64(reinterpret_cast<std::size_t>(complex_num.img));
            return dd::combineHash(h1, h2);
        }
    };
} // namespace std

#endif //DD_PACKAGE_COMPLEX_HPP
