/*
 * Copyright (c) 2017, K. Kumar, Delft University of Technology (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <functional>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <astro/astro.hpp>
#include <integrate/integrate.hpp>

#include "dustsim/dynamicalSystem.hpp"
#include "dustsim/outputWriter.hpp"
#include "dustsim/singleParticleSimulator.hpp"
#include "dustsim/tools.hpp"

namespace dustsim
{

//! Execute single particle simulator.
void executeSingleParticleSimulator( const rapidjson::Document& config )
{
    // Verify config parameters. Exception is thrown if any of the parameters are missing.
    const SingleParticleSimulatorInput input = checkSingleParticleSimulatorInput( config );

    std::cout << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << "                           Run simulator                          " << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << std::endl;

    // Compute initial state in Cartesian elements.
    const State initialState = astro::convertKeplerianToCartesianElements(
        input.initialStateKeplerianElements, input.gravitationalParameter );
    std::cout << "Cartesian initial state            (" << initialState << ")" << std::endl;
    std::cout << std::endl;

    // Set current tune and state to initial time and state.
    State currentState = initialState;
    Real currentTime = input.startTime;

    // Create instance of dynamical system.
    std::cout << "Setting up dynamical model ..." << std::endl;
    DynamicalSystem dynamics( input.gravitationalParameter,
                              input.isJ2AccelerationModelActive,
                              input.j2Coefficient,
                              input.equatorialRadius,
                              input.isRadiationPressureAccelerationModelActive,
                              input.particleRadius,
                              input.particleBulkDensity,
                              input.radiationPressureCoefficient );
    std::cout << "Dynamical model set up successfully!" << std::endl;
    std::cout << std::endl;

    // Write metadata to file.
    std::ofstream metadataFile( input.metadataFilePath );
    print( metadataFile, "gravitational_parameter", input.gravitationalParameter, "km^{3} s^{-2}" );
    metadataFile << std::endl;
    print( metadataFile, "j2_coefficient", input.j2Coefficient, "-" );
    metadataFile << std::endl;
    print( metadataFile, "equatorial_radius", input.equatorialRadius, "km" );
    metadataFile << std::endl;
    std::ostringstream initialKeplerStateString;
    initialKeplerStateString << input.initialStateKeplerianElements[ astro::semiMajorAxisIndex ];
    initialKeplerStateString << "; ";
    initialKeplerStateString << input.initialStateKeplerianElements[ astro::eccentricityIndex ];
    initialKeplerStateString << "; ";
    initialKeplerStateString
        << sml::convertRadiansToDegrees(
            input.initialStateKeplerianElements[ astro::inclinationIndex ] );
    initialKeplerStateString << "; ";
    initialKeplerStateString
        << sml::convertRadiansToDegrees(
            input.initialStateKeplerianElements[ astro::argumentOfPeriapsisIndex ] );
    initialKeplerStateString << "; ";
    initialKeplerStateString
        << sml::convertRadiansToDegrees(
            input.initialStateKeplerianElements[ astro::longitudeOfAscendingNodeIndex ] );
    initialKeplerStateString << "; ";
    initialKeplerStateString
        << sml::convertRadiansToDegrees(
            input.initialStateKeplerianElements[ astro::trueAnomalyIndex ] );
    print( metadataFile, "initial_state_kepler", initialKeplerStateString.str( ), "km | deg" );
    metadataFile << std::endl;
    print( metadataFile, "start_time", input.startTime, "s" );
    metadataFile << std::endl;
    print( metadataFile, "end_time", input.endTime, "s" );
    metadataFile << std::endl;
    print( metadataFile, "time_step", input.timeStep, "s" );
    metadataFile << std::endl;

    // Create file stream to write state history to.
    std::ofstream stateHistoryFile( input.stateHistoryFilePath );
    stateHistoryFile << "t,x,y,z,xdot,ydot,zdot,a,e,i,aop,raan,ta" << std::endl;
    StateHistoryWriter writer( stateHistoryFile,
                               input.gravitationalParameter );

    // Set up numerical integrator.
    std::cout << "Exeucting numerical integration ..." << std::endl;
    if ( input.integrator == rk4 )
    {
        auto stateDerivativePointer = std::bind( &DynamicalSystem::operator( ),
                                                 &dynamics,
                                                 std::placeholders::_1,
                                                 std::placeholders::_2 );

        while ( currentTime < input.endTime )
        {
            integrate::stepRK4< Real, State >( currentTime,
                                               currentState,
                                               input.timeStep,
                                               stateDerivativePointer );
            writer( currentTime, currentState );
        }
    }

    std::cout << "Numerical integrator executed successfully!" << std::endl;
}

//! Check input parameters for single_particle_simulator application mode.
SingleParticleSimulatorInput checkSingleParticleSimulatorInput( const rapidjson::Document& config )
{
    // Extract central gravity model parameters.
    const Real gravitationalParameter
        = find( config, "gravitational_parameter" )->value.GetDouble( );
    std::cout << "Gravitational parameter            " << gravitationalParameter
              << " [km^3 s^-2]" << std::endl;

    // Extract J2 gravity model parameters.
    const bool j2AcclerationModelFlag = find( config, "is_j2_active" )->value.GetBool( );
    std::cout << "Is J2 model active?                " << j2AcclerationModelFlag << std::endl;

    const Real j2Coefficient = find( config, "j2_coefficient" )->value.GetDouble( );
    std::cout << "J2 coefficient                     " << j2Coefficient << " [-]" << std::endl;

    const Real equatorialRadius
        = find( config, "equatorial_radius" )->value.GetDouble( );
    std::cout << "Equatorial radius                  " << equatorialRadius << " [km]" << std::endl;

    // Extract radiation pressure model parameters.
    const bool radiationPressureFlag
        = find( config, "is_radiation_pressure_active" )->value.GetBool( );
    std::cout << "Is SRP model active?               "
              << ( radiationPressureFlag ? "true" : "false" ) << std::endl;

    const Real particleRadius = find( config, "particle_radius" )->value.GetDouble( );
    std::cout << "Particle radius                    "
              << particleRadius << " [micron]" << std::endl;

    const Real particleBulkDensity = find( config, "particle_bulk_density" )->value.GetDouble( );
    std::cout << "Particle bulk density              "
              << particleBulkDensity << " [kg m^-3]" << std::endl;

    const Real radiationPressureCoefficient
        = find( config, "radiation_pressure_coefficient" )->value.GetDouble( );
    std::cout << "Radiation pressure coefficient     "
              << radiationPressureCoefficient << " [-]" << std::endl;

    // Extract initial state of dust particle in Keplerian elements.
    ConfigIterator initialStateKeplerianElementsIterator = find( config, "initial_state_kepler" );
    Vector initialStateKeplerianElementsVector( 6 );

    initialStateKeplerianElementsVector[ astro::semiMajorAxisIndex ]
     = initialStateKeplerianElementsIterator->value[ astro::semiMajorAxisIndex ].GetDouble( );
    initialStateKeplerianElementsVector[ astro::eccentricityIndex ]
     = initialStateKeplerianElementsIterator->value[ astro::eccentricityIndex ].GetDouble( );
    initialStateKeplerianElementsVector[ astro::inclinationIndex ]
     = sml::convertDegreesToRadians(
         initialStateKeplerianElementsIterator->value[ astro::inclinationIndex ].GetDouble( ) );
    initialStateKeplerianElementsVector[ astro::argumentOfPeriapsisIndex ]
     = sml::convertDegreesToRadians(
         initialStateKeplerianElementsIterator->value[
             astro::argumentOfPeriapsisIndex ].GetDouble( ) );
    initialStateKeplerianElementsVector[ astro::longitudeOfAscendingNodeIndex ]
     = sml::convertDegreesToRadians(
         initialStateKeplerianElementsIterator->value[
             astro::longitudeOfAscendingNodeIndex ].GetDouble( ) );
    initialStateKeplerianElementsVector[ astro::trueAnomalyIndex ]
     = sml::convertDegreesToRadians(
         initialStateKeplerianElementsIterator->value[ astro::trueAnomalyIndex ].GetDouble( ) );

    const State initialStateKeplerianElements( initialStateKeplerianElementsVector );

    std::cout << "Initial state (Kepler)             ("
           << initialStateKeplerianElements[ astro::semiMajorAxisIndex ] << ", "
           << initialStateKeplerianElements[ astro::eccentricityIndex ] << ", "
           << initialStateKeplerianElements[ astro::inclinationIndex ] << ", "
           << initialStateKeplerianElements[ astro::argumentOfPeriapsisIndex ] << ", "
           << initialStateKeplerianElements[ astro::longitudeOfAscendingNodeIndex ] << ", "
           << initialStateKeplerianElements[ astro::trueAnomalyIndex ] << ") "
           << "[km, -, rad, rad, rad, rad]" << std::endl;

    // Extract selected numerical integrator.
    const std::string integratorString = find( config, "integrator" )->value.GetString( );
    Integrator integrator = rk4;
    if ( integratorString.compare( "rk4" ) != 0 )
    {
        if ( integratorString.compare( "rkf78" ) == 0 )
        {
            integrator = rkf78;
        }
        else if ( integratorString.compare( "dopri5" ) == 0 )
        {
            integrator = dopri5;
        }
        else if ( integratorString.compare( "bs" ) == 0 )
        {
            integrator = bs;
        }
        else
        {
            std::cout << std::endl;
            std::cerr << "Selected numerical integrator \""
                      << integratorString
                      << "\" is incorrect!" << std::endl;
            throw;
        }
    }
    std::cout << "Integrator                         " << integratorString << std::endl;

    // Extract integrator time settings.
    const Real stateTime            = find( config, "start_time" )->value.GetDouble( );
    std::cout << "Start epoch                        " << stateTime << " [s]" << std::endl;
    const Real endTime              = find( config, "end_time" )->value.GetDouble( );
    std::cout << "End epoch                          " << endTime << " [s]" << std::endl;
    const Real timeStep             = find( config, "time_step" )->value.GetDouble( );
    std::cout << "Time step                          " << timeStep << " [s]" << std::endl;

    // Extract integrator tolerances.
    const Real relativeTolerance    = find( config, "relative_tolerance" )->value.GetDouble( );
    std::cout << "Relative tolerance                 " << relativeTolerance << " [-]" << std::endl;
    const Real absoluteTolerance    = find( config, "absolute_tolerance" )->value.GetDouble( );
    std::cout << "Absolute tolerance                 " << absoluteTolerance << " [-]" << std::endl;

    // Extract file writer settings.
    const std::string metadataFilePath
        = find( config, "metadata_file_path" )->value.GetString( );
    std::cout << "Metadata file path                 " << metadataFilePath << std::endl;
    const std::string stateHistoryFilePath
        = find( config, "state_history_file_path" )->value.GetString( );
    std::cout << "State history file path            " << stateHistoryFilePath << std::endl;

    return SingleParticleSimulatorInput( gravitationalParameter,
                                         j2AcclerationModelFlag,
                                         j2Coefficient,
                                         equatorialRadius,
                                         radiationPressureFlag,
                                         particleRadius,
                                         particleBulkDensity,
                                         radiationPressureCoefficient,
                                         initialStateKeplerianElements,
                                         integrator,
                                         stateTime,
                                         endTime,
                                         timeStep,
                                         relativeTolerance,
                                         absoluteTolerance,
                                         metadataFilePath,
                                         stateHistoryFilePath );
}

} // namespace dustsim
