/*
 * Copyright (c) 2017, K. Kumar, Delft University of Technology (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <fstream>
#include <iostream>

#include <boost/numeric/odeint.hpp>

#include <astro/astro.hpp>

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
    State initialState = astro::convertKeplerianToCartesianElements(
        input.initialStateKeplerianElements, input.gravitationalParameter );
    std::cout << "Cartesian initial state       (";
    for ( unsigned int i = 0; i < initialState.size( ) - 1; i++ )
    {
        std::cout << initialState[ i ] << ", ";
    }
    std::cout << initialState[ initialState.size( ) - 1 ] << ")" << std::endl;
    std::cout << std::endl;

    // Set current state to initial state.
    State currentState = initialState;

    // // Set current epoch to start epoch.
    // Real currentEpoch = input.startEpoch;

    // Create instance of dynamical system.
    std::cout << "Setting up dynamical model ..." << std::endl;
    DynamicalSystem dynamics( input.gravitationalParameter,
                              input.j2Coefficient,
                              input.equatorialRadius );
    std::cout << "Dynamical model set up successfully!" << std::endl;
    std::cout << std::endl;

    // Create file stream to write state history to.
    std::ofstream stateHistoryFile( input.stateHistoryFilePath );
    stateHistoryFile << "t,x,y,z,xdot,ydot,zdot,a,e,i,aop,raan,ta" << std::endl;
    StateHistoryWriter writer( stateHistoryFile,
                               input.gravitationalParameter );

    // Set up numerical integrator.
    std::cout << "Setting up numerical integrator ..." << std::endl;
    boost::numeric::odeint::integrate_const( boost::numeric::odeint::runge_kutta4< State >(),
                                             dynamics,
                                             currentState,
                                             input.startEpoch,
                                             input.endEpoch,
                                             input.timeStep,
                                             writer  );
    std::cout << "Numerical integrator set up successfully!" << std::endl;
}

//! Check input parameters for single_particle_simulator application mode.
SingleParticleSimulatorInput checkSingleParticleSimulatorInput( const rapidjson::Document& config )
{
    // Extract environment model parameters.
    const Real gravitationalParameter
        = find( config, "gravitational_parameter" )->value.GetDouble( );
    std::cout << "Gravitational parameter       " << gravitationalParameter
              << " [km^3 s^-2]" << std::endl;

    const Real j2Coefficient = find( config, "j2_coefficient" )->value.GetDouble( );
    std::cout << "J2 coefficient                " << j2Coefficient << " [-]" << std::endl;

    const Real equatorialRadius
        = find( config, "equatorial_radius" )->value.GetDouble( );
    std::cout << "Equatorial radius             " << equatorialRadius << " [km]" << std::endl;

   // Extract initial state of dust particle in Keplerian elements.
    ConfigIterator initialStateKeplerianElementsIterator = find( config, "initial_state_kepler" );
    State initialStateKeplerianElements;
    initialStateKeplerianElements[ astro::semiMajorAxisIndex ]
     = initialStateKeplerianElementsIterator->value[ astro::semiMajorAxisIndex ].GetDouble( );
    initialStateKeplerianElements[ astro::eccentricityIndex ]
     = initialStateKeplerianElementsIterator->value[ astro::eccentricityIndex ].GetDouble( );
    initialStateKeplerianElements[ astro::inclinationIndex ]
     = sml::convertDegreesToRadians(
         initialStateKeplerianElementsIterator->value[ astro::inclinationIndex ].GetDouble( ) );
    initialStateKeplerianElements[ astro::argumentOfPeriapsisIndex ]
     = sml::convertDegreesToRadians(
         initialStateKeplerianElementsIterator->value[
             astro::argumentOfPeriapsisIndex ].GetDouble( ) );
    initialStateKeplerianElements[ astro::longitudeOfAscendingNodeIndex ]
     = sml::convertDegreesToRadians(
         initialStateKeplerianElementsIterator->value[
             astro::longitudeOfAscendingNodeIndex ].GetDouble( ) );
    initialStateKeplerianElements[ astro::trueAnomalyIndex ]
     = sml::convertDegreesToRadians(
         initialStateKeplerianElementsIterator->value[ astro::trueAnomalyIndex ].GetDouble( ) );
    std::cout << "Initial state (Kepler)        ("
           << initialStateKeplerianElements[ astro::semiMajorAxisIndex ] << ", "
           << initialStateKeplerianElements[ astro::eccentricityIndex ] << ", "
           << initialStateKeplerianElements[ astro::inclinationIndex ] << ", "
           << initialStateKeplerianElements[ astro::argumentOfPeriapsisIndex ] << ", "
           << initialStateKeplerianElements[ astro::longitudeOfAscendingNodeIndex ] << ", "
           << initialStateKeplerianElements[ astro::trueAnomalyIndex ] << ") "
           << "[km, -, rad, rad, rad, rad]" << std::endl;

     // Extract integrator time settings.
    const Real startEpoch           = find( config, "start_epoch" )->value.GetDouble( );
    std::cout << "Start epoch                   " << startEpoch << " [s]" << std::endl;
    const Real endEpoch             = find( config, "end_epoch" )->value.GetDouble( );
    std::cout << "End epoch                     " << endEpoch << " [s]" << std::endl;
    const Real timeStep             = find( config, "time_step" )->value.GetDouble( );
    std::cout << "Time step                     " << timeStep << " [s]" << std::endl;

    // Extract integrator tolerances.
    const Real relativeTolerance    = find( config, "relative_tolerance" )->value.GetDouble( );
    std::cout << "Relative tolerance            " << relativeTolerance << " [-]" << std::endl;
    const Real absoluteTolerance    = find( config, "absolute_tolerance" )->value.GetDouble( );
    std::cout << "Absolute tolerance            " << absoluteTolerance << " [-]" << std::endl;

    // Extract file writer settings.
    const std::string stateHistoryFilePath
        = find( config, "state_history_file_path" )->value.GetString( );
    std::cout << "State history file path       " << stateHistoryFilePath << std::endl;

    return SingleParticleSimulatorInput( gravitationalParameter,
                                         j2Coefficient,
                                         equatorialRadius,
                                         initialStateKeplerianElements,
                                         startEpoch,
                                         endEpoch,
                                         timeStep,
                                         relativeTolerance,
                                         absoluteTolerance,
                                         stateHistoryFilePath );
}

} // namespace dustsim
