# Copyright (c) 2017, K. Kumar, Delft University of Technology (me@kartikkumar.com)
# Distributed under the MIT License.
# See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT

# Include script to build external libraries with CMake.
include(ExternalProject)

# -------------------------------

# Boost: https://boost.org

find_package(Boost)

if(NOT APPLE)
  include_directories(SYSTEM AFTER "${Boost_INCLUDE_DIRS}")
else(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${Boost_INCLUDE_DIRS}\"")
endif(NOT APPLE)

# -------------------------------

# SQLiteCpp: https://github.com/srombauts/sqlitecpp

if(NOT BUILD_DEPENDENCIES)
  find_package(SQLiteCpp 2002000)
endif(NOT BUILD_DEPENDENCIES)

if(NOT SQLITECPP_FOUND)
  message(STATUS "SQLiteCpp will be downloaded when ${CMAKE_PROJECT_NAME} is built")
  ExternalProject_Add(sqlitecpp-lib
    PREFIX ${EXTERNAL_PATH}/SQLiteCpp
    #--Download step--------------
    URL https://github.com/SRombauts/SQLiteCpp/archive/master.zip
    TIMEOUT 30
    #--Update/Patch step----------
    #--Configure step-------------
    #--Build step-----------------
    BUILD_IN_SOURCE 1
    #--Install step---------------
    INSTALL_COMMAND ""
    #--Output logging-------------
    LOG_DOWNLOAD ON
  )
  ExternalProject_Get_Property(sqlitecpp-lib source_dir)
  set(SQLITECPP_INCLUDE_DIRS ${source_dir}/include
    CACHE INTERNAL "Path to include folder for SQLiteCpp")
  set(SQLITECPP_LIBRARY_DIR ${source_dir} CACHE INTERNAL "Path to library folder for SQLiteCpp")
  set(SQLITECPP_LIBRARY "SQLiteCpp")
endif(NOT SQLITECPP_FOUND)

if(NOT APPLE)
  include_directories(SYSTEM AFTER "${SQLITECPP_INCLUDE_DIRS}")
else(NOT APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SQLITECPP_INCLUDE_DIRS}\"")
endif(NOT APPLE)
link_directories(${SQLITECPP_LIBRARY_DIR})

# -------------------------------

# rapidjson: https://github.com/miloyip/rapidjson

if(NOT BUILD_DEPENDENCIES)
  find_package(rapidjson)
endif(NOT BUILD_DEPENDENCIES)

if(NOT RAPIDJSON_FOUND)
  message(STATUS "RapidJSON will be downloaded when ${CMAKE_PROJECT_NAME} is built")
  ExternalProject_Add(rapidjson-lib
    PREFIX ${EXTERNAL_PATH}/RapidJson
    #--Download step--------------
    URL https://github.com/miloyip/rapidjson/archive/master.zip
    TIMEOUT 30
    #--Update/Patch step----------
    #--Configure step-------------
    CONFIGURE_COMMAND ""
    #--Build step-----------------
    BUILD_COMMAND ""
    #--Install step---------------
    INSTALL_COMMAND ""
    #--Output logging-------------
    LOG_DOWNLOAD ON
  )
  ExternalProject_Get_Property(rapidjson-lib source_dir)
  set(RAPIDJSON_INCLUDE_DIRS ${source_dir}/include
    CACHE INTERNAL "Path to include folder for RapidJSON")
endif(NOT RAPIDJSON_FOUND)

if(NOT APPLE)
  include_directories(SYSTEM AFTER "${RAPIDJSON_INCLUDE_DIRS}")
else(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${RAPIDJSON_INCLUDE_DIRS}\"")
endif(NOT APPLE)

# -------------------------------

# sml: https://github.com/openastro/sml

if(NOT BUILD_DEPENDENCIES)
  find_package(sml)
endif(NOT BUILD_DEPENDENCIES)

if(NOT SML_FOUND)
  message(STATUS "sml will be downloaded when ${CMAKE_PROJECT_NAME} is built")
  ExternalProject_Add(sml-lib
    PREFIX ${EXTERNAL_PATH}/sml
    #--Download step--------------
    URL https://github.com/openastro/sml/archive/master.zip
    TIMEOUT 30
    #--Update/Patch step----------
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    #--Configure step-------------
    CONFIGURE_COMMAND ""
    #--Build step-----------------
    BUILD_COMMAND ""
    #--Install step---------------
    INSTALL_COMMAND ""
    #--Output logging-------------
    LOG_DOWNLOAD ON
  )
  ExternalProject_Get_Property(sml-lib source_dir)
  set(SML_INCLUDE_DIRS ${source_dir}/include CACHE INTERNAL "Path to include folder for sml")
endif(NOT SML_FOUND)

if(NOT APPLE)
  include_directories(SYSTEM AFTER "${SML_INCLUDE_DIRS}")
else(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${SML_INCLUDE_DIRS}\"")
endif(NOT APPLE)

# -------------------------------

# astro: https://github.com/openastro/astro

if(NOT BUILD_DEPENDENCIES)
  find_package(astro)
endif(NOT BUILD_DEPENDENCIES)

if(NOT ASTRO_FOUND)
  message(STATUS "astro will be downloaded when ${CMAKE_PROJECT_NAME} is built")
  ExternalProject_Add(astro-lib
    DEPENDS sml-lib
    PREFIX ${EXTERNAL_PATH}/astro
    #--Download step--------------
    URL https://github.com/openastro/astro/archive/master.zip
    TIMEOUT 30
    #--Update/Patch step----------
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    #--Configure step-------------
    CONFIGURE_COMMAND ""
    #--Build step-----------------
    BUILD_COMMAND ""
    #--Install step---------------
    INSTALL_COMMAND ""
    #--Output logging-------------
    LOG_DOWNLOAD ON
  )
  ExternalProject_Get_Property(astro-lib source_dir)
  set(ASTRO_INCLUDE_DIRS ${source_dir}/include CACHE INTERNAL "Path to include folder for astro")
endif(NOT ASTRO_FOUND)

if(NOT APPLE)
  include_directories(SYSTEM AFTER "${ASTRO_INCLUDE_DIRS}")
else(APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${ASTRO_INCLUDE_DIRS}\"")
endif(NOT APPLE)

# -------------------------------

# Catch: https://github.com/philsquared/Catch

if(BUILD_TESTS)
  if(NOT BUILD_DEPENDENCIES)
    find_package(CATCH)
  endif(NOT BUILD_DEPENDENCIES)

  if(NOT CATCH_FOUND)
    message(STATUS "Catch will be downloaded when ${CMAKE_PROJECT_NAME} is built")
    ExternalProject_Add(catch-lib
      PREFIX ${EXTERNAL_PATH}/Catch
      #--Download step--------------
      URL https://github.com/philsquared/Catch/archive/master.zip
      TIMEOUT 30
      #--Update/Patch step----------
      UPDATE_COMMAND ""
      PATCH_COMMAND ""
      #--Configure step-------------
      CONFIGURE_COMMAND ""
      #--Build step-----------------
      BUILD_COMMAND ""
      #--Install step---------------
      INSTALL_COMMAND ""
      #--Output logging-------------
      LOG_DOWNLOAD ON
    )
    ExternalProject_Get_Property(catch-lib source_dir)
    set(CATCH_INCLUDE_DIRS ${source_dir}/include CACHE INTERNAL "Path to include folder for Catch")
  endif(NOT CATCH_FOUND)

  if(NOT APPLE)
    include_directories(SYSTEM AFTER "${CATCH_INCLUDE_DIRS}")
  else(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"${CATCH_INCLUDE_DIRS}\"")
  endif(NOT APPLE)
endif(BUILD_TESTS)

# -------------------------------
