/*
 * Copyright (c) 2009-2022 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include <sml/smlAll.hpp>
#include <astro/astroAll.hpp>
#include <integrate/integrateAll.hpp>

#include "dustsim/dynamicalSystem.hpp"
#include "dustsim/singleParticleSimulator.hpp"
#include "dustsim/tools.hpp"

namespace dustsim
{

//! Execute single particle simulator.
void executeSingleParticleSimulator(const nlohmann::json& config)
{
    // Verify config parameters. Exception is thrown if any of the parameters are missing.
    const SingleParticleSimulatorInput input = checkSingleParticleSimulatorInput(config);

    std::cout << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << "                 Compute additional model parameters              " << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << std::endl;

    // Compute mean motion of central body around the Sun [rad/s].
    const Real solarMeanMotion = astro::computeKeplerMeanMotion(
        input.solarDistance * astro::ASTRO_AU_IN_KM, input.solarGravitationalParameter);
    std::cout << "  Solar mean motion                " << solarMeanMotion
              << " [rad s^-1]" << std::endl;

    // Compute radiation pressure for complete absorption at distance of central body from the
    // Sun [N m^-2].
    const Real radiationPressure = astro::computeAbsorptionRadiationPressure(input.solarEnergyFlux);
    std::cout << "  Radiation pressure               " << radiationPressure
              << " [N m^2]" << std::endl;

    // Compute initial state in Cartesian elements.
    const State initialState = astro::convertKeplerianToCartesianElements(
        input.initialStateKeplerianElements, input.gravitationalParameter);
    std::cout << "  Cartesian initial state          (" << initialState << ")" << std::endl;

    std::cout << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << "                           Run simulator                          " << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << std::endl;

    // Set current time and state and step size to initial values specified by user.
    std::cout << "Setting up numerical integrator settings ..." << std::endl;
    State state = initialState;
    Real time = input.startTime;
    Real stepSize = input.timeStep;
    std::cout << "Numerical integrator settings set up successfully!" << std::endl;
    std::cout << std::endl;

    // Create instance of dynamical system.
    std::cout << "Setting up dynamical model ..." << std::endl;
    DynamicalSystem dynamics(input.gravitationalParameter,
                             input.isJ2AccelerationModelActive,
                             input.j2Coefficient,
                             input.equatorialRadius,
                             input.isRadiationPressureAccelerationModelActive,
                             input.particleRadius,
                             input.particleBulkDensity,
                             input.radiationPressureCoefficient,
                             solarMeanMotion,
                             radiationPressure);
    std::cout << "Dynamical model set up successfully!" << std::endl;
    std::cout << std::endl;

    // Write metadata to file.
    std::ofstream metadataFile(input.metadataFilePath);
    print(metadataFile, "gravitational_parameter", input.gravitationalParameter, "km^{3} s^{-2}");
    metadataFile << std::endl;
    print(metadataFile, "j2_coefficient", input.j2Coefficient, "-");
    metadataFile << std::endl;
    print(metadataFile, "equatorial_radius", input.equatorialRadius, "km");
    metadataFile << std::endl;
    std::ostringstream initialKeplerStateString;
    initialKeplerStateString << input.initialStateKeplerianElements[astro::semiMajorAxisIndex];
    initialKeplerStateString << "; ";
    initialKeplerStateString << input.initialStateKeplerianElements[astro::eccentricityIndex];
    initialKeplerStateString << "; ";
    initialKeplerStateString
        << sml::convertRadiansToDegrees(
            input.initialStateKeplerianElements[astro::inclinationIndex]);
    initialKeplerStateString << "; ";
    initialKeplerStateString
        << sml::convertRadiansToDegrees(
            input.initialStateKeplerianElements[astro::argumentOfPeriapsisIndex]);
    initialKeplerStateString << "; ";
    initialKeplerStateString
        << sml::convertRadiansToDegrees(
            input.initialStateKeplerianElements[astro::longitudeOfAscendingNodeIndex]);
    initialKeplerStateString << "; ";
    initialKeplerStateString
        << sml::convertRadiansToDegrees(
            input.initialStateKeplerianElements[astro::trueAnomalyIndex]);
    print(metadataFile, "initial_state_kepler", initialKeplerStateString.str(), "km | deg");
    metadataFile << std::endl;
    print(metadataFile, "start_time", input.startTime, "s");
    metadataFile << std::endl;
    print(metadataFile, "end_time", input.endTime, "s");
    metadataFile << std::endl;
    print(metadataFile, "time_step", input.timeStep, "s");
    metadataFile << std::endl;
    print(metadataFile, "solar_mean_motion", solarMeanMotion, "rad s^{-1}");
    metadataFile << std::endl;
    print(metadataFile, "radiation_pressure", radiationPressure, "N m^{2}");
    metadataFile << std::endl;

    // Create file stream to write state history to.
    std::ofstream stateHistoryFile(input.stateHistoryFilePath);
    stateHistoryFile << "t,dt,x,y,z,xdot,ydot,zdot,a,e,i,aop,raan,ta" << std::endl;

    // Write initial data to state history file.
    stateHistoryFile << std::setprecision(std::numeric_limits<double>::digits10)
                     << time << ","
                     << stepSize << ", "
                     << state << ","
                     << input.initialStateKeplerianElements << std::endl;

    // Set up member function pointer to point to state derivative function defined for
    // the dynamical system.
    auto stateDerivativePointer = std::bind(&DynamicalSystem::operator(),
                                            &dynamics,
                                            std::placeholders::_1,
                                            std::placeholders::_2);

    // Set up numerical integrator.
    std::cout << "Executing numerical integration and writing results to file ..." << std::endl;
    if (input.integrator == rk4)
    {
        while (time < input.endTime)
        {
            integrate::stepRK4<Real, State>(time,
                                            state,
                                            stepSize,
                                            stateDerivativePointer);

            const State stateInKeplerElements
                = astro::convertCartesianToKeplerianElements(state,
                                                             input.gravitationalParameter);

            stateHistoryFile << std::setprecision(std::numeric_limits<double>::digits10)
                             << time << ","
                             << stepSize << ","
                             << state << ","
                             << stateInKeplerElements << std::endl;
        }
    }
    // // // else if (input.integrator == rkf45)
    // // // {
    // // //     while (time < input.endTime)
    // // //     {
    // // //         integrate::stepRKF45<Real, State>(time,
    // // //                                           state,
    // // //                                           stepSize,
    // // //                                           stateDerivativePointer,
    // // //                                           input.relativeTolerance,
    // // //                                           input.minimumStepSize,
    // // //                                           input.maximumStepSize);

    // // //         const State stateInKeplerElements
    // // //             = astro::convertCartesianToKeplerianElements(state,
    // // //                                                          input.gravitationalParameter);

    // // //         stateHistoryFile << std::setprecision(std::numeric_limits<double>::digits10)
    // // //                          << time << ","
    // // //                          << stepSize << ","
    // // //                          << state << ","
    // // //                          << stateInKeplerElements << std::endl;
        // }
    // }
    else if (input.integrator == rkf78)
    {
    // // //     Int outputIntervalCounter = 1;

    // // //     while (time < input.endTime)
    // // //     {
    // // //         Real previousTime = time;
    // // //         State previousState = state;

    // // //         integrate::stepRKF78<Real, State>(time,
    // // //                                           state,
    // // //                                           stepSize,
    // // //                                           stateDerivativePointer,
    // // //                                           input.relativeTolerance,
    // // //                                           input.minimumStepSize,
    // // //                                           input.maximumStepSize);

    // // //         Real nextOutputTime = input.startTime + input.outputInterval * outputIntervalCounter;

    // // //         if (input.outputInterval > 0.0)
    // // //         {
    // // //             while (time > nextOutputTime + input.minimumStepSize)
    // // //             {
    // // //                 Real outputStepSize = (nextOutputTime - previousTime);

    // // //                 integrate::stepRKF78<Real, State>(previousTime,
    // // //                                                   previousState,
    // // //                                                   outputStepSize,
    // // //                                                   stateDerivativePointer,
    // // //                                                   input.relativeTolerance,
    // // //                                                   input.minimumStepSize,
    // // //                                                   input.maximumStepSize);
    // // //                 const State stateInKeplerElements
    // // //                     = astro::convertCartesianToKeplerianElements(previousState,
    // // //                                                                   input.gravitationalParameter);

    // // //                 stateHistoryFile << std::setprecision(std::numeric_limits<double>::digits10)
    // // //                                  << previousTime << ","
    // // //                                  << input.outputInterval << ","
    // // //                                  << previousState << ","
    // // //                                  << stateInKeplerElements << std::endl;

    // // //                 outputIntervalCounter++;
    // // //                 nextOutputTime = input.startTime + input.outputInterval * outputIntervalCounter;
    // // //             }
    // // //         }
    // // //         else
    // // //         {
    // // //             const State stateInKeplerElements
    // // //                 = astro::convertCartesianToKeplerianElements(state,
    // // //                                                               input.gravitationalParameter);

    // // //             stateHistoryFile << std::setprecision(std::numeric_limits<double>::digits10)
    // // //                              << time << ","
    // // //                              << stepSize << ","
    // // //                              << state << ","
    // // //                              << stateInKeplerElements << std::endl;
    // // //         }
    // // //     }
    }

    std::cout << "Numerical integrator executed successfully and results written to file!"
              << std::endl;
}

//! Check input parameters for single_particle_simulator application mode.
SingleParticleSimulatorInput checkSingleParticleSimulatorInput(const nlohmann::json& config)
{
    // Extract central gravity model parameters.
    const Real gravitationalParameter = config.at("gravitational_parameter").get<Real>();
    std::cout << "  Gravitational parameter          " << gravitationalParameter
              << " [km^3 s^-2]" << std::endl;

    // Extract J2 gravity model parameters.
    const bool j2AcclerationModelFlag = config.at("is_j2_active").get<bool>();
    std::cout << "  Is J2 model active?              " << j2AcclerationModelFlag << std::endl;

    const Real j2Coefficient = config.at("j2_coefficient").get<Real>();
    std::cout << "  J2 coefficient                   " << j2Coefficient << " [-]" << std::endl;

    const Real equatorialRadius = config.at("equatorial_radius").get<Real>();
    std::cout << "  Equatorial radius                " << equatorialRadius << " [km]" << std::endl;

    // Extract radiation pressure model parameters.
    const bool radiationPressureFlag
        = config.at("is_radiation_pressure_active").get<bool>();
    std::cout << "  Is SRP model active?             "
              << (radiationPressureFlag ? "true" : "false") << std::endl;

    const Real particleRadius = config.at("particle_radius").get<Real>() * 1.0e-6;
    std::cout << "  Particle radius                  " << particleRadius << " [m]" << std::endl;

    const Real particleBulkDensity = config.at("particle_bulk_density").get<Real>();
    std::cout << "  Particle bulk density            "
              << particleBulkDensity << " [kg m^-3]" << std::endl;

    const Real radiationPressureCoefficient
        = config.at("radiation_pressure_coefficient").get<Real>();
    std::cout << "  Radiation pressure coefficient   "
              << radiationPressureCoefficient << " [-]" << std::endl;

    const Real solarDistance = config.at("mean_solar_distance").get<Real>();
    std::cout << "  Mean solar distance              " << solarDistance << " [AU]" << std::endl;

    const Real solarGravitationalParameter
        = config.at("solar_gravitatational_parameter").get<Real>();
    std::cout << "  Solar gravitational parameter    "
              << solarGravitationalParameter << " [km^3 s^-2]" << std::endl;

    const Real solarEnergyFlux  = config.at("mean_solar_energy_flux").get<Real>();
    std::cout << "  Mean solar energy flux           "
              << solarEnergyFlux << " [W m^-2]" << std::endl;

    // Extract initial state of dust particle in Keplerian elements.
    // ConfigIterator initialStateKeplerianElementsIterator = config.at("initial_state_kepler");
    Vector initialStateKeplerianElementsVector(6);

    initialStateKeplerianElementsVector[astro::semiMajorAxisIndex]
        = config.at("initial_state_kepler")[astro::semiMajorAxisIndex].get<Real>();
    initialStateKeplerianElementsVector[astro::eccentricityIndex]
        = config.at("initial_state_kepler")[astro::eccentricityIndex].get<Real>();
    initialStateKeplerianElementsVector[astro::inclinationIndex]
     = sml::convertDegreesToRadians(
         config.at("initial_state_kepler")[astro::inclinationIndex].get<Real>());
    initialStateKeplerianElementsVector[astro::argumentOfPeriapsisIndex]
     = sml::convertDegreesToRadians(
         config.at("initial_state_kepler")[
             astro::argumentOfPeriapsisIndex].get<Real>());
    initialStateKeplerianElementsVector[astro::longitudeOfAscendingNodeIndex]
     = sml::convertDegreesToRadians(
         config.at("initial_state_kepler")[
             astro::longitudeOfAscendingNodeIndex].get<Real>());
    initialStateKeplerianElementsVector[astro::trueAnomalyIndex]
     = sml::convertDegreesToRadians(
         config.at("initial_state_kepler")[astro::trueAnomalyIndex].get<Real>());

    const State initialStateKeplerianElements(initialStateKeplerianElementsVector);

    std::cout << "  Initial state (Kepler)           ("
              << initialStateKeplerianElements[astro::semiMajorAxisIndex] << ", "
              << initialStateKeplerianElements[astro::eccentricityIndex] << ", "
              << initialStateKeplerianElements[astro::inclinationIndex] << ", "
              << initialStateKeplerianElements[astro::argumentOfPeriapsisIndex] << ", "
              << initialStateKeplerianElements[astro::longitudeOfAscendingNodeIndex] << ", "
              << initialStateKeplerianElements[astro::trueAnomalyIndex] << ") "
              << "[km, -, rad, rad, rad, rad]" << std::endl;

    // // Extract selected numerical integrator.
    const std::string integratorString = config.at("integrator").get<std::string>();
    Integrator integrator = rk4;
    if (integratorString.compare("rk4") != 0)
    {
        if (integratorString.compare("rkf45") == 0)
        {
            integrator = rkf45;
        }
        else if (integratorString.compare("rkf78") == 0)
        {
            integrator = rkf78;
        }
        else
        {
            std::cout << std::endl;
            std::cerr << "Selected numerical integrator \""
                      << integratorString
                      << "\" is unavailable!" << std::endl;
            throw;
        }
    }
    std::cout << "  Integrator                       " << integratorString << std::endl;

    // Extract integrator time settings.
    const Real stateTime            = config.at("start_time").get<Real>();
    std::cout << "  Start epoch                      " << stateTime << " [s]" << std::endl;
    const Real endTime              = config.at("end_time").get<Real>();
    std::cout << "  End epoch                        " << endTime << " [s]" << std::endl;
    const Real timeStep             = config.at("time_step").get<Real>();
    std::cout << "  Time step                        " << timeStep << " [s]" << std::endl;

    // Extract integrator tolerances.
    const Real relativeTolerance    = config.at("relative_tolerance").get<Real>();
    std::cout << "  Relative tolerance               " << relativeTolerance << " [-]" << std::endl;
    const Real absoluteTolerance    = config.at("absolute_tolerance").get<Real>();
    std::cout << "  Absolute tolerance               " << absoluteTolerance << " [-]" << std::endl;

    // Extract integrator step size bounds.
    const Real minimumStepSize    = config.at("minimum_step_size").get<Real>();
    std::cout << "  Minimum step size                " << minimumStepSize << " [s]" << std::endl;
    const Real maximumStepSize    = config.at("maximum_step_size").get<Real>();
    std::cout << "  Maximum step size                " << maximumStepSize << " [s]" << std::endl;

    const Real outputInterval    = config.at("output_interval").get<Real>();
    std::cout << "  Output interval                  " << outputInterval << " [s]" << std::endl;

    // Extract file writer settings.
    const std::string rootDirectory = config.at("root_directory").get<std::string>();
    const std::string metadataFilePath
        = rootDirectory + "/" + config.at("metadata_file_path").get<std::string>();
    std::cout << "  Metadata file path               " << metadataFilePath << std::endl;
    const std::string stateHistoryFilePath
        = rootDirectory + "/" + config.at("state_history_file_path").get<std::string>();
    std::cout << "  State history file path          " << stateHistoryFilePath << std::endl;

    return SingleParticleSimulatorInput(gravitationalParameter,
                                        j2AcclerationModelFlag,
                                        j2Coefficient,
                                        equatorialRadius,
                                        radiationPressureFlag,
                                        particleRadius,
                                        particleBulkDensity,
                                        radiationPressureCoefficient,
                                        solarDistance,
                                        solarGravitationalParameter,
                                        solarEnergyFlux,
                                        initialStateKeplerianElements,
                                        integrator,
                                        stateTime,
                                        endTime,
                                        timeStep,
                                        relativeTolerance,
                                        absoluteTolerance,
                                        minimumStepSize,
                                        maximumStepSize,
                                        outputInterval,
                                        metadataFilePath,
                                        stateHistoryFilePath);
}

} // namespace dustsim
