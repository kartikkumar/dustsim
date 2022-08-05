# Copyright (c) 2009-2022 Kartik Kumar (me@kartikkumar.com)
# Distributed under the MIT License.
# See accompanying file LICENSE or copy at http://opensource.org/licenses/MIT

# The CMake setup for this project is based off of the following source:
# - https://cliutils.gitlab.io/modern-cmake

# Make sure that an old version of CMake is not being used so modern features are available
# Update this to the mix, max versions that are applicable
cmake_minimum_required(VERSION 3.14)
if(${CMAKE_VERSION} VERSION_LESS 3.23)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
    cmake_policy(VERSION 3.23)
endif()

# Set up project details
project(
    dustsim
    VERSION     2.0
    DESCRIPTION "A circumplanetary dust simulator"
    LANGUAGES   CXX
    )

# Require out-of-source builds
# Source: https://stackoverflow.com/questions/11143062/getting-cmake-to-build-out-of-source-without-wrapping-scripts
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt.")
endif()

# Support folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Enable downloads during the configure step using FetchContent
# FetchContent was added in CMake 3.11
# FetchContent_MakeAvailable was not added until CMake 3.14
include(FetchContent)

# Fetch JSON library
# Source: https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
FetchContent_Declare(json
  GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
  GIT_PROGRESS TRUE
  GIT_SHALLOW TRUE
  GIT_TAG v3.10.5)
FetchContent_MakeAvailable(json)

# Fetch sml library
FetchContent_Declare(sml
  GIT_REPOSITORY https://github.com/openastro/sml
  GIT_PROGRESS TRUE
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(sml)

# Fetch integrate library
FetchContent_Declare(integrate
  GIT_REPOSITORY https://github.com/openastro/integrate
  GIT_PROGRESS TRUE
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(integrate)

# Fetch astro library
FetchContent_Declare(astro
  GIT_REPOSITORY https://github.com/openastro/astro
  GIT_PROGRESS TRUE
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(astro)

# Find SQLite3 library
find_package(SQLite3)

# Fetch SQLiteCpp library
FetchContent_Declare(SQLiteCpp
  GIT_REPOSITORY https://github.com/SRombauts/SQLiteCpp
  GIT_PROGRESS TRUE
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(SQLiteCpp)

# Instantiate version file based on project details
configure_file (
    "${PROJECT_SOURCE_DIR}/version.hpp.in"
    "${PROJECT_BINARY_DIR}/include/${project_name}/version.hpp"
)

# Add subdirs for headers, sources, and executables
# Subdirs contain CMakeLists.txt files with commands to setup build
# Comment the following line if the library is header-only
add_subdirectory(src)
# Uncomment the following line if the library is header-only
#add_subdirectory(include)
# Comment the following line if there are no applications in this project
add_subdirectory(apps)

# Enable testing
if(BUILD_TESTING)
    include(CTest)
    add_subdirectory(tests)
endif(BUILD_TESTING)

# Build Doxygen docs
if(BUILD_DOCS)
    find_package(Doxygen OPTIONAL_COMPONENTS dot)
    if(Doxygen_FOUND)
        add_subdirectory(docs)
    else()
        message(STATUS "Doxygen not found, not building docs")
    endif(Doxygen_FOUND)
endif(BUILD_DOCS)

# Set up uninstall target to enable 'make uninstall'
# Source: https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
