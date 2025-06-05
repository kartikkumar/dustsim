/*
 * Copyright (c) 2009-2025 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

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
    // Extract root directory where input & output files are stored
    // N.B.: directory has to already exist!
    std::string ioDirectory = config.at("io_directory").get<std::string>();
    std::cout << "I/O directory                    " << ioDirectory << std::endl;

    // Verify config parameters. Exception is thrown if any of the parameters are missing.
    const SingleParticleSimulatorInput input = checkSingleParticleSimulatorInput(config);

    std::cout << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << "                 Compute additional model parameters              " << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << std::endl;

    // Compute initial state in Cartesian elements.
    const State initialState = astro::convertKeplerianToCartesianElements(
        input.initialStateKeplerianElements, input.gravitationalParameter);
    std::cout << "Cartesian initial state          (" << initialState << ")"
              << "[km, km, km, km/s, kms, km/s]" << std::endl;

    // Compute mean motion of central body around the Sun [rad/s].
    const Real solarMeanMotion = astro::computeKeplerMeanMotion(
        input.solarDistance * astro::ASTRO_AU_IN_KM, input.solarGravitationalParameter);
    if (input.isRadiationPressureAccelerationModelActive)
    {
        std::cout << "Solar mean motion                " << solarMeanMotion
                  << " [rad s^-1]" << std::endl;
    }

    // Compute radiation pressure for complete absorption at distance of central body from the
    // Sun [N m^-2].
    const Real radiationPressure = astro::computeAbsorptionRadiationPressure(input.solarEnergyFlux);

    if (input.isRadiationPressureAccelerationModelActive)
    {
        std::cout << "Radiation pressure               " << radiationPressure
                  << " [N m^2]" << std::endl;
    }

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

    if (input.isJ2AccelerationModelActive)
    {
        print(metadataFile, "is_j2_active", 'y', "");
        metadataFile << std::endl;
        print(metadataFile, "equatorial_radius", input.equatorialRadius, "km");
        metadataFile << std::endl;
        print(metadataFile, "j2_coefficient", input.j2Coefficient, "-");
        metadataFile << std::endl;
    }
    else
    {
        print(metadataFile, "is_j2_active", 'n', "");
        metadataFile << std::endl;
        print(metadataFile, "equatorial_radius", '-', "km");
        metadataFile << std::endl;
        print(metadataFile, "j2_coefficient", '-', "-");
        metadataFile << std::endl;
    }

    if (input.isRadiationPressureAccelerationModelActive)
    {
        print(metadataFile, "is_radiation_pressure_active", 'y', "");
        metadataFile << std::endl;
        print(metadataFile, "particle_radius", input.particleRadius, "10^{-6} m");
        metadataFile << std::endl;
        print(metadataFile, "particle_bulk_density", input.particleBulkDensity, "kg m^{-3}");
        metadataFile << std::endl;
        print(metadataFile, "radiation_pressure_coefficient",
            input.radiationPressureCoefficient, "kg m^{-3}");
        metadataFile << std::endl;
        print(metadataFile, "mean_solar_distance", input.solarDistance, "rad s^{-1}");
        metadataFile << std::endl;
        print(metadataFile, "solar_gravitatational_parameter",
            input.solarGravitationalParameter, "N m^{2}");
        metadataFile << std::endl;
        print(metadataFile, "mean_solar_energy_flux", input.solarEnergyFlux, "N m^{2}");
        metadataFile << std::endl;
    }
    else
    {
        print(metadataFile, "is_radiation_pressure_active", 'n', "");
        metadataFile << std::endl;
        print(metadataFile, "particle_radius", "-", "10^{-6} m");
        metadataFile << std::endl;
        print(metadataFile, "particle_bulk_density", "-", "kg m^{-3}");
        metadataFile << std::endl;
        print(metadataFile, "radiation_pressure_coefficient", "-", "kg m^{-3}");
        metadataFile << std::endl;
        print(metadataFile, "mean_solar_distance", "-", "rad s^{-1}");
        metadataFile << std::endl;
        print(metadataFile, "solar_gravitatational_parameter", "-", "N m^{2}");
        metadataFile << std::endl;
        print(metadataFile, "mean_solar_energy_flux", "-", "N m^{2}");
        metadataFile << std::endl;
    }

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
    if (input.integrator == rk4)
    {
        print(metadataFile, "integrator", "rk4", "");
    }
    else if (input.integrator == rkf45)
    {
        print(metadataFile, "integrator", "rkf45", "");
    }
    else if (input.integrator == rkf78)
    {
        print(metadataFile, "integrator", "rkf78", "");
    }
    else
    {
        print(metadataFile, "integrator", "error", "");
    }
    metadataFile << std::endl;
    print(metadataFile, "start_time", input.startTime, "s");
    metadataFile << std::endl;
    print(metadataFile, "end_time", input.endTime, "s");
    metadataFile << std::endl;
    print(metadataFile, "time_step", input.timeStep, "s");
    metadataFile << std::endl;
    if (input.integrator == rkf45)
    {
        print(metadataFile, "relative_tolerance", input.relativeTolerance, "-");
        metadataFile << std::endl;
        print(metadataFile, "absolute_tolerance", input.absoluteTolerance, "-");
        metadataFile << std::endl;
        print(metadataFile, "minimum_step_size", input.minimumStepSize, "s");
        metadataFile << std::endl;
        print(metadataFile, "maximum_step_size", input.maximumStepSize, "s");
        metadataFile << std::endl;
    }

    // Create file stream to write state history to.
    std::ofstream stateHistoryFile(input.stateHistoryFilePath);
    stateHistoryFile
        << "t_s,dt_s,x_km,y_km,z_km,xdot_km_s,ydot_km_s,"
        << "zdot_km_s,a_km,e,i_rad,aop_rad,raan_rad,ta_rad"
        << std::endl;

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
    else if (input.integrator == rkf45)
    {
        while (time < input.endTime)
        {
            integrate::stepRKF45<Real, State>(time,
                                              state,
                                              stepSize,
                                              stateDerivativePointer,
                                              input.relativeTolerance,
                                              input.minimumStepSize,
                                              input.maximumStepSize);

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
    else if (input.integrator == rkf78)
    {
        // Int outputIntervalCounter = 1;

        while (time < input.endTime)
        {
            // Real previousTime = time;
            // State previousState = state;

            integrate::stepRKF78<Real, State>(time,
                                              state,
                                              stepSize,
                                              stateDerivativePointer,
                                              input.relativeTolerance,
                                              input.minimumStepSize,
                                              input.maximumStepSize);

            // Real nextOutputTime = input.startTime + input.outputInterval * outputIntervalCounter;

            if (config.contains("output_interval"))
            {
                // while (time > nextOutputTime + input.minimumStepSize)
                // {
                //     Real outputStepSize = (nextOutputTime - previousTime);

                //     integrate::stepRKF78<Real, State>(previousTime,
                //                                       previousState,
                //                                       outputStepSize,
                //                                       stateDerivativePointer,
                //                                       input.relativeTolerance,
                //                                       input.minimumStepSize,
                //                                       input.maximumStepSize);
                //     const State stateInKeplerElements
                //         = astro::convertCartesianToKeplerianElements(previousState,
                //                                                       input.gravitationalParameter);

                //     stateHistoryFile << std::setprecision(std::numeric_limits<double>::digits10)
                //                      << previousTime << ","
                //                      << input.outputInterval << ","
                //                      << previousState << ","
                //                      << stateInKeplerElements << std::endl;

                //     outputIntervalCounter++;
                //     nextOutputTime = input.startTime + input.outputInterval * outputIntervalCounter;
                // }
            }
            else
            {
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
    }

    std::cout << "Numerical integrator executed successfully and results written to file!"
              << std::endl;
}

//! Check input parameters for single particle simulator application mode.
SingleParticleSimulatorInput checkSingleParticleSimulatorInput(const nlohmann::json& config)
{
    // Extract central gravity model parameters.
    const Real gravitationalParameter = config["gravitational_parameter"].get<Real>();
    std::cout << "Gravitational parameter          " << gravitationalParameter
              << " [km^3 s^-2]" << std::endl;

    const Real equatorialRadius = config["equatorial_radius"].get<Real>();
    std::cout << "Equatorial radius                " << equatorialRadius << " [km]" << std::endl;

    // Extract J2 gravity model parameters.
    const bool j2AcclerationModelFlag = config["is_j2_active"].get<bool>();
    std::cout << "Is J2 model active?              "
              << (j2AcclerationModelFlag ? "true" : "false") << std::endl;

    const Real j2Coefficient = config["j2_coefficient"].get<Real>();
    if (j2AcclerationModelFlag)
    {
        std::cout << "J2 coefficient                   " << j2Coefficient << " [-]" << std::endl;
    }

    // Extract radiation pressure model parameters.
    const bool radiationPressureFlag
        = config["is_radiation_pressure_active"].get<bool>();
    std::cout << "Is SRP model active?             "
              << (radiationPressureFlag ? "true" : "false") << std::endl;

    const Real particleRadius = config["particle_radius"].get<Real>() * 1.0e-6;
    const Real particleBulkDensity = config["particle_bulk_density"].get<Real>();
    const Real radiationPressureCoefficient
        = config["radiation_pressure_coefficient"].get<Real>();
    const Real solarDistance = config["mean_solar_distance"].get<Real>();
    const Real solarGravitationalParameter
        = config["solar_gravitatational_parameter"].get<Real>();
    const Real solarEnergyFlux  = config["mean_solar_energy_flux"].get<Real>();

    if (radiationPressureFlag)
    {
        std::cout << "Particle radius                  " << particleRadius << " [m]" << std::endl;

        std::cout << "Particle bulk density            "
                  << particleBulkDensity << " [kg m^-3]" << std::endl;
        std::cout << "Radiation pressure coefficient   "
                  << radiationPressureCoefficient << " [-]" << std::endl;
        std::cout << "Mean solar distance              " << solarDistance << " [AU]" << std::endl;
        std::cout << "Solar gravitational parameter    "
                  << solarGravitationalParameter << " [km^3 s^-2]" << std::endl;
        std::cout << "Mean solar energy flux           "
                  << solarEnergyFlux << " [W m^-2]" << std::endl;
    }

    // Extract initial state of dust particle in Keplerian elements.
    const State initialStateKeplerianElements(
      { config["initial_state_kepler"][astro::semiMajorAxisIndex].get<Real>(),
        config["initial_state_kepler"][astro::eccentricityIndex].get<Real>(),
        sml::convertDegreesToRadians(
          config["initial_state_kepler"][
            astro::inclinationIndex].get<Real>()),
        sml::convertDegreesToRadians(
          config["initial_state_kepler"][
            astro::argumentOfPeriapsisIndex].get<Real>()),
        sml::convertDegreesToRadians(
          config["initial_state_kepler"][
            astro::longitudeOfAscendingNodeIndex].get<Real>()),
        sml::convertDegreesToRadians(
          config["initial_state_kepler"][
            astro::trueAnomalyIndex].get<Real>())});

    std::cout << "Initial state (Kepler)           ("
              << initialStateKeplerianElements[astro::semiMajorAxisIndex] << ", "
              << initialStateKeplerianElements[astro::eccentricityIndex] << ", "
              << config["initial_state_kepler"][astro::inclinationIndex].get<Real>() << ", "
              << config["initial_state_kepler"][
                astro::argumentOfPeriapsisIndex].get<Real>() << ", "
              << config["initial_state_kepler"][
                astro::longitudeOfAscendingNodeIndex].get<Real>() << ", "
              << config["initial_state_kepler"][astro::trueAnomalyIndex].get<Real>() << ") "
              << "[km, -, deg, deg, deg, deg]" << std::endl;

    // Extract selected numerical integrator.
    const std::string integratorString = config["integrator"].get<std::string>();
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
    std::cout << "Integrator                       " << integratorString << std::endl;

    // Extract integrator time settings.
    const Real stateTime            = config["start_time"].get<Real>();
    std::cout << "Start epoch                      " << stateTime << " [s]" << std::endl;
    const Real endTime              = config["end_time"].get<Real>();
    std::cout << "End epoch                        " << endTime << " [s]" << std::endl;
    const Real timeStep             = config["time_step"].get<Real>();
    std::cout << "Time step                        " << timeStep << " [s]" << std::endl;

    // Extract integrator tolerances & step size bounds.
    const Real relativeTolerance    = config["relative_tolerance"].get<Real>();
    const Real absoluteTolerance    = config["absolute_tolerance"].get<Real>();
    const Real minimumStepSize      = config["minimum_step_size"].get<Real>();
    const Real maximumStepSize      = config["maximum_step_size"].get<Real>();
    Real outputInterval = 0.0;
    if (config.contains("output_interval"))
    {
        outputInterval              = config["output_interval"].get<Real>();
    }
    else
    {
        outputInterval              = std::numeric_limits<Real>::quiet_NaN();
    }

    if (integrator != rk4)
    {
        std::cout << "Relative tolerance               "
            << relativeTolerance << " [-]" << std::endl;
        std::cout << "Absolute tolerance               "
            << absoluteTolerance << " [-]" << std::endl;
        std::cout << "Minimum step size                " << minimumStepSize << " [s]" << std::endl;
        std::cout << "Maximum step size                " << maximumStepSize << " [s]" << std::endl;
        std::cout << "Output interval                  " << outputInterval << " [s]" << std::endl;
    }

    // Extract file writer settings.
    const std::string ioDirectory = config["io_directory"].get<std::string>();
    const std::string metadataFilePath
        = ioDirectory + "/" + config["metadata_file_path"].get<std::string>();
    std::cout << "Metadata file path               " << metadataFilePath << std::endl;
    const std::string stateHistoryFilePath
        = ioDirectory + "/" + config["state_history_file_path"].get<std::string>();
    std::cout << "State history file path          " << stateHistoryFilePath << std::endl;

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
