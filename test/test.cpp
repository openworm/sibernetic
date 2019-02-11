#include "solver_container.hpp"
#include <gtest/gtest.h>

int add(int i, int j) { return i + j; }

TEST(test, test) { EXPECT_EQ(2, add(1, 1)); }

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}