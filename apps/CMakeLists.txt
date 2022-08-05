# Copyright (c) 2009-2022 Kartik Kumar (me@kartikkumar.com)
# Distributed under the MIT License.
# See accompanying file LICENSE or copy at http://opensource.org/licenses/MIT

# The CMake setup for this project is based off of the following source:
# - https://cliutils.gitlab.io/modern-cmake

# Create executables and link necessary libraries for each application

# Set up main application
add_executable(dustsim dustsim.cpp)
target_compile_features(dustsim PRIVATE cxx_std_11)
target_link_libraries(dustsim PRIVATE dustsim_lib nlohmann_json::nlohmann_json)
