#include "dd/MDDPackage.hpp"

#include "gtest/gtest.h"

#include <cstddef>

TEST(DDPackageTest, TrivialTest) {
  auto dd = std::make_unique<dd::MDDPackage>(2, std::vector<std::size_t>{2,3});
  EXPECT_EQ(dd->qregisters(), 2);
}
