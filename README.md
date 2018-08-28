dustsim
===

\cond [![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT) \endcond

dustsim is a circumplanetary dust simulator that can be used to study dust dynamics, as a result of gravitational and non-gravitational forces acting on a dust particle.

Features
------

  - General directory structure common to C++ projects
  - Install script (`make install`)
  - CPack script for packaging (`make package`)
  - Automatic documentation generation ([Doxygen](http://www.doxygen.org "Doxygen homepage"))
  - Separate file to specify location of project files (`ProjectFiles.cmake`)

Requirements
------

To install this project, please ensure that you have installed the following (install guides are provided on the respective websites):

  - [Git](http://git-scm.com)
  - A C++ compiler, e.g., [GCC](https://gcc.gnu.org/), [clang](http://clang.llvm.org/), [MinGW](http://www.mingw.org/)
  - [CMake](http://www.cmake.org "CMake homepage")
  - [Doxygen](http://www.doxygen.org "Doxygen homepage") (optional)

Installation
------

Run the following commands to download, build, and install this project. The ` --depth 1` parameter passed to `git clone` ensures that the git history is not downloaded. In case you would like to preserve the history of this project, omit that option.

    git clone https://www.github.com/kartikkumar/dustsim --depth 1
    cd dustsim
    git submodule init && git submodule update
    mkdir build && cd build
    cmake .. && cmake --build .

To push this project to your own remote repository, you can run the following command, which will overcome the issues with utilizing a shallow clone:

    git commit --amend .

This rewrites the last commit and ensures that you can then push the repository to a remote (e.g., Github, BitBucket, Gitlab, etc.).

To install the header files, run the following from within the `build` directory:

    make install

Note that dependencies are installed by fetching them online, in case they cannot be detected on your local system. If the build process fails, check the error log given. Typically, building fails due to timeout. Simply run the `cmake --build .` command once more.

Build options
-------------

You can pass the following, general command-line options when running CMake:

  - `-DCMAKE_INSTALL_PREFIX[=$install_dir]`: set path prefix for install script (`make install`); if not set, defaults to usual locations
  - `-DBUILD_SHARED_LIBS=[ON|OFF (default)]`: build shared libraries instead of static
  - `-DBUILD_DOXYGEN_DOCS[=ON|OFF (default)]`: build the [Doxygen](http://www.doxygen.org "Doxygen homepage") documentation ([LaTeX](http://www.latex-project.org/) must be installed with `amsmath` package)
  - `-DBUILD_DEPENDENCIES[=ON|OFF (default)]`: force local build of dependencies, instead of first searching system-wide using `find_package()`

Project structure
-------------

This project has been set up with a specific file/folder structure in mind. The following describes some important features of this setup:

  - `cmake/Modules`: Contains `CMake` modules
  - `config`: Contains templates for configuration files to run dustsim (for C++ simulator and Python plotting scripts)
  - `docs`: Contains code documentation generated by [Doxygen](http://www.doxygen.org "Doxygen homepage")
  - `include`: Project header files (*.hpp)
  - `python`: Python plotting scripts for data generated by simulator
  - `src`: Project source files (*.cpp), including `main.cpp`
  - `CMakeLists.txt`: main `CMakelists.txt` file for project (should not need to be modified for basic build)
  - `Dependencies.cmake`: list of dependencies and automated build, triggered if dependency cannot be found locally
  - `Doxyfile.in`: [Doxygen](http://www.doxygen.org "Doxygen homepage") configuration file
  - `LICENSE.md`: license file for project
  - `ProjectFiles.cmake`: list of project source files to build
  - `README.md`: project readme file

Contributing
------------

Once you've made your great commits:

  1. [Fork](https://github.com/kartikkumar/dustsim/fork) dustsim
  2. Create a topic branch - `git checkout -b my_branch`
  3. Push to your branch - `git push origin my_branch`
  4. Create a [Pull Request](http://help.github.com/pull-requests/) from your branch
  5. That's it!

Disclaimer
------------

The copyright holders are not liable for any damage(s) incurred due to improper use of dustsim.

TODO
------------

@todo Find a way to have nested variables in `Doxygen` config file so that e.g., `@@CMAKE_PROJECT_NAME@_VERSION@` works.

@todo Add version detection in `CMake` module so that find_package respects minimum version required.

@todo Find a way to provide an option to clean installation.

@todo Find a way to remove \cond \endcond workaround to get Doxygen to not throw warnings in readme.
