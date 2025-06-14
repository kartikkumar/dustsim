/*
 * Copyright (c) 2009-2025 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <nlohmann/json.hpp>

#include "dustsim/singleParticleSimulator.hpp"
// #include "dustsim/bulkParticleSimulator.hpp"

int main(const int numberOfInputs, const char* inputArguments[])
{
	///////////////////////////////////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "                              dustsim                             " << std::endl;
    std::cout << std::endl;
    std::cout << "      Copyright (c) 2009-2025, K. Kumar (me@kartikkumar.com)      " << std::endl;
    std::cout << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;

    /////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << "                       Parse input parameters                     " << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << std::endl;

    // Check that only one input has been provided (a JSON file).
    if (numberOfInputs - 1 != 1)
    {
        std::cerr << "ERROR: Number of inputs is wrong. Please only provide a JSON input file!"
                  << std::endl;
        throw;
    }

    // Read and store JSON input document (filter out comment lines).
    // TODO: Need to make comment-line filtering more robust.
    std::ifstream inputFile(inputArguments[1]);
    std::stringstream jsonDocumentBuffer;
    std::string inputLine;
    while (std::getline(inputFile, inputLine))
    {
        size_t startPosition = inputLine.find_first_not_of(" \t");
        if (std::string::npos != startPosition)
        {
            inputLine = inputLine.substr(startPosition);
        }

        if (inputLine.substr(0, 2) != "//")
        {
            jsonDocumentBuffer << inputLine << "\n";
        }
    }

    nlohmann::json config = nlohmann::json::parse(jsonDocumentBuffer.str().c_str());

    if (!config.contains("mode"))
    {
        std::cerr << "ERROR: Configuration option \"mode\" could not be found in JSON input!"
                  << std::endl;
        throw;
    }

    std::string mode = config.at("mode").get<std::string>();
    if (mode.compare("single_particle_simulator") == 0)
    {
        std::cout << "Mode                             " << mode << std::endl;
        dustsim::executeSingleParticleSimulator(config);
    }
    else if (mode.compare("bulk_particle_simulator") == 0)
    {
        std::cout << "  Mode                             " << mode << std::endl;
    // //     // dustsim::executeBulkParticleSimulator(config);
    }
    else
    {
        std::cout << std::endl;
        std::cerr << "ERROR: Requested mode \"" << mode << "\" is invalid!" << std::endl;
        throw;
    }

    /////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    std::cout << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "                        Exited successfully!                      " << std::endl;
    std::cout << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;

    return EXIT_SUCCESS;

    ///////////////////////////////////////////////////////////////////////////
}
