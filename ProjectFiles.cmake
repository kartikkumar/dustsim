# Copyright (c) 2009-2018, K. Kumar (me@kartikkumar.com)
# Distributed under the MIT License.
# See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT

# Set project source files.
set(SRC
  "${SRC_PATH}/bulkParticleSimulator.cpp"
  "${SRC_PATH}/singleParticleSimulator.cpp"
  "${SRC_PATH}/state.cpp"
  "${SRC_PATH}/tools.cpp"
)

# Set project main file.
set(MAIN_SRC
  "${SRC_PATH}/main.cpp"
)

# Set project test source files.
set(TEST_SRC
  "${TEST_SRC_PATH}/testDustsim.cpp"
  "${TEST_SRC_PATH}/testDummy.cpp"
)
