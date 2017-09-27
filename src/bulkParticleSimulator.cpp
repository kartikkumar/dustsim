/*
 * Copyright (c) 2017, K. Kumar, Delft University of Technology (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <SQLiteCpp/SQLiteCpp.h>

#include "dustsim/bulkParticleSimulator.hpp"
#include "dustsim/tools.hpp"

namespace dustsim
{

//! Execute bulk particle simulator.
void executeBulkParticleSimulator( const rapidjson::Document& config )
{
    // Verify config parameters. Exception is thrown if any of the parameters are missing.
    const BulkParticleSimulatorInput input = checkBulkParticleSimulatorInput( config );

    // Open database in read/write mode.
    SQLite::Database database( input.databaseFilePath, SQLite::OPEN_READWRITE|SQLite::OPEN_CREATE );

    // Initialize random generator.
    boost::random::mt19937 randomGenerator;

    // Set up generators for initial states for dust particles.
    boost::random::uniform_real_distribution< Real > semiMajorAxisGenerator(
        input.semiMajorAxisMinimum, input.semiMajorAxisMaximum );

    for ( int i = 0; i < input.numberOfParticles; ++i )
    {
        std::cout << std::setprecision( 15 ) << semiMajorAxisGenerator( randomGenerator ) << std::endl;
    }
}

//! Check input parameters for bulk_particle_simulator application mode.
BulkParticleSimulatorInput checkBulkParticleSimulatorInput( const rapidjson::Document& config )
{
    // Extract central gravity model parameters.
    const Real gravitationalParameter
        = find( config, "gravitational_parameter" )->value.GetDouble( );
    std::cout << "Gravitational parameter            " << gravitationalParameter
              << " [km^3 s^-2]" << std::endl;

    // Extract J2 gravity model parameters.
    const bool j2AcclerationModelFlag = find( config, "is_j2_active" )->value.GetBool( );
    std::cout << "Is J2 model active?                "
              << ( j2AcclerationModelFlag ? "true" : "false" ) << std::endl;

    const Real j2Coefficient = find( config, "j2_coefficient" )->value.GetDouble( );
    std::cout << "J2 coefficient                     " << j2Coefficient << " [-]" << std::endl;

    const Real equatorialRadius
        = find( config, "equatorial_radius" )->value.GetDouble( );
    std::cout << "Equatorial radius                  " << equatorialRadius << " [km]" << std::endl;

    // Extract solar radiation pressure model parameters.
    const bool solarRadiationPressureFlag
        = find( config, "is_solar_radiation_pressure_active" )->value.GetBool( );
    std::cout << "Is SRP model active?               "
              << ( solarRadiationPressureFlag ? "true" : "false" ) << std::endl;

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

    // Extract number of particles to simulate.
    const Int numberOfParticles = find( config, "number_of_particles" )->value.GetInt( );
    std::cout << "Number of particles                " << numberOfParticles << std::endl;

    // Extract parameters for distribution of initial states of dust particle in Keplerian elements.
    const Real semiMajorAxisMinimum = find( config, "semi_major_axis_minimum" )->value.GetDouble( );
    std::cout << "Semi-major axis minimum            "
              << semiMajorAxisMinimum << " [km]" << std::endl;

    const Real semiMajorAxisMaximum = find( config, "semi_major_axis_maximum" )->value.GetDouble( );
    std::cout << "Semi-major axis maximum            "
              << semiMajorAxisMaximum << " [km]" << std::endl;

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
    const Real startEpoch           = find( config, "start_epoch" )->value.GetDouble( );
    std::cout << "Start epoch                        " << startEpoch << " [s]" << std::endl;
    const Real endEpoch             = find( config, "end_epoch" )->value.GetDouble( );
    std::cout << "End epoch                          " << endEpoch << " [s]" << std::endl;
    const Real timeStep             = find( config, "time_step" )->value.GetDouble( );
    std::cout << "Time step                          " << timeStep << " [s]" << std::endl;

    // Extract integrator tolerances.
    const Real relativeTolerance    = find( config, "relative_tolerance" )->value.GetDouble( );
    std::cout << "Relative tolerance                 " << relativeTolerance << " [-]" << std::endl;
    const Real absoluteTolerance    = find( config, "absolute_tolerance" )->value.GetDouble( );
    std::cout << "Absolute tolerance                 " << absoluteTolerance << " [-]" << std::endl;

    // Extract SQLite database settings.
    const std::string databaseFilePath = find( config, "database_file_path" )->value.GetString( );
    std::cout << "Database file path                 " << databaseFilePath << std::endl;

    return BulkParticleSimulatorInput( gravitationalParameter,
                                       j2AcclerationModelFlag,
                                       j2Coefficient,
                                       equatorialRadius,
                                       numberOfParticles,
                                       semiMajorAxisMinimum,
                                       semiMajorAxisMaximum,
                                       integrator,
                                       startEpoch,
                                       endEpoch,
                                       timeStep,
                                       relativeTolerance,
                                       absoluteTolerance,
                                       databaseFilePath );
}

} // namespace dustsim