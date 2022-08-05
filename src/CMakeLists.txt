# Copyright (c) 2009-2022 Kartik Kumar (me@kartikkumar.com)
# Distributed under the MIT License.
# See accompanying file LICENSE or copy at http://opensource.org/licenses/MIT

# The CMake setup for this project is based off of the following source:
# - https://cliutils.gitlab.io/modern-cmake

# -----------------------------------------------

# List all files that should be included in the library here
set(
  SOURCE_LIST
  singleParticleSimulator.cpp
  state.cpp
  )

# -----------------------------------------------

# Grab a list of all header files
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/dustsim/*.hpp")

# Make an automatic library - will be static or dynamic based on user setting
# Note that headers are optional, and do not affect add_library, but they will not
# show up in IDEs unless they are listed in add_library
add_library(
  dustsim_lib
  ${SOURCE_LIST}
  ${HEADER_LIST}
)
target_include_directories(dustsim_lib PUBLIC ../include)
target_compile_features(dustsim_lib PRIVATE cxx_std_11)
target_link_libraries(dustsim_lib PRIVATE sml_lib integrate_lib astro_lib SQLiteCpp nlohmann_json::nlohmann_json)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST}
  )
