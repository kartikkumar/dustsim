/*
 * Copyright (c) 2017, K. Kumar, Delft University of Technology (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

#include <boost/random/uniform_real_distribution.hpp>

#include <keplerian_toolbox.h>

#include <SQLiteCpp/SQLiteCpp.h>

#include <astro/astro.hpp>

#include <sml/sml.hpp>

#include "dustsim/bulkParticleSimulator.hpp"
#include "dustsim/tools.hpp"

namespace dustsim
{

//! Execute bulk particle simulator.
void executeBulkParticleSimulator( const rapidjson::Document& config )
{
    // Verify config parameters. Exception is thrown if any of the parameters are missing.
    const BulkParticleSimulatorInput input = checkBulkParticleSimulatorInput( config );

    std::cout << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << "                           Run simulator                          " << std::endl;
    std::cout << "******************************************************************" << std::endl;
    std::cout << std::endl;

    std::cout << "Setting up random number generators for initial states ..." << std::endl;

    // Initialize random seed.
    boost::random::mt19937 randomSeed;

    // Set up generators for initial states for dust particles.
    // Semi-major axis values are generated from a uniform distribution, with minimum and maximum
    // set in the config file.
    boost::random::uniform_real_distribution< Real > semiMajorAxisGenerator(
        input.semiMajorAxisMinimum, input.semiMajorAxisMaximum );

    // Eccentricity & argument of periapses values are generated from normal distributions for the
    // components of the eccentricity vector [h_e = e*cos( AoP ), k_e = e*sin( AoP )].
    const Real eccentricityStandardDeviation
        = input.eccentricityFullWidthHalfMaximum / ( 2.0 * std::sqrt( 2.0 * std::log( 2 ) ) );
    std::cout << "Eccentricity standard deviation    " << eccentricityStandardDeviation
              << " [-]" << std::endl;

    boost::random::normal_distribution< Real > eccentricityComponentGenerator(
        0.0, eccentricityStandardDeviation );

    // Inclination & longitude of ascending node values are generated from normal distributions for
    // the components of the inclination vector [h_i = i*cos( RAAN ), k_i = i*sin( RAAN ) ].
    const Real inclinationStandardDeviation
        = input.inclinationFullWidthHalfMaximum / ( 2.0 * std::sqrt( 2.0 * std::log( 2 ) ) );
    std::cout << "Inclination standard deviation     " << inclinationStandardDeviation
              << " [rad]" << std::endl;

    boost::random::normal_distribution< Real > inclinationComponentGenerator(
        0.0, inclinationStandardDeviation );

    // Mean anomaly values are generated from a uniform distribution.
    boost::random::uniform_real_distribution< Real > meanAnomalyGenerator( 0.0, 2.0 * sml::SML_PI );

    std::cout << "Random number generators set up successfully!" << std::endl;
    std::cout << std::endl;

    std::cout << "Setting up database ..." << std::endl;
    std::cout << std::endl;

    // Open database in read/write mode.
    SQLite::Database database( input.databaseFilePath, SQLite::OPEN_READWRITE|SQLite::OPEN_CREATE );

    // Set names of tables in databaes.
    std::string initialStatesTableName = "initial_states";

    std::cout << "Creating table '" << initialStatesTableName << "' ..." << std::endl;

    std::ostringstream initialStatesDropTable;
    initialStatesDropTable << "DROP TABLE IF EXISTS " << initialStatesTableName << ";";
    database.exec( initialStatesDropTable.str( ) );

    // Create input states table.
    std::ostringstream initialStatesTableCreate;
    initialStatesTableCreate
        << "CREATE TABLE initial_states ("
        << "\"simulation_id\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
        << "\"simulated\" INTEGER NOT NULL,"
        << "\"semi_major_axis\" REAL NOT NULL,"
        << "\"eccentricity\" REAL NOT NULL,"
        << "\"inclination\" REAL NOT NULL,"
        << "\"argument_of_periapsis\" REAL NOT NULL,"
        << "\"longitude_of_ascending_node\" REAL NOT NULL,"
        << "\"true_anomaly\" REAL NOT NULL);";

    database.exec( initialStatesTableCreate.str( ) );

    std::cout << "Database set up successfully!" << std::endl;
    std::cout << std::endl;

    std::cout << "Populating initial states table ..." << std::endl;

    std::ostringstream initialStatesTableInsert;
    initialStatesTableInsert
        << "INSERT INTO '" << initialStatesTableName << "' VALUES ("
        << "NULL,"
        << ":simulated,"
        << ":semi_major_axis,"
        << ":eccentricity,"
        << ":inclination,"
        << ":argument_of_periapsis,"
        << ":longitude_of_ascending_node,"
        << ":true_anomaly"
        << ");";

    SQLite::Statement query( database, initialStatesTableInsert.str( ) );

    // Set up database transaction.
    SQLite::Transaction initialStatesTableInsertTransaction( database );

    for ( int i = 0; i < input.numberOfParticles ; ++i )
    {
        // Generate random semi-major axis [km].
        const Real semiMajorAxis = semiMajorAxisGenerator( randomSeed );

        // Generate random eccentricity vector components.
        const Real eccentricityXComponent = eccentricityComponentGenerator( randomSeed );
        const Real eccentricityYComponent = eccentricityComponentGenerator( randomSeed );

        // Compute eccentricity [-] and argument of periapsis [rad].
        const Real eccentricity = std::sqrt( eccentricityXComponent * eccentricityXComponent
                                             + eccentricityYComponent * eccentricityYComponent );

        const Real argumentOfPeriapsis
            = std::atan2( eccentricityYComponent, eccentricityXComponent );

        // Generate random inclination vector components.
        const Real inclinationXComponent = inclinationComponentGenerator( randomSeed );
        const Real inclinationYComponent = inclinationComponentGenerator( randomSeed );

        // Compute inclination [rad] and longitude of ascending node [rad].
        const Real inclination = std::sqrt( inclinationXComponent * inclinationXComponent
                                             + inclinationYComponent * inclinationYComponent );

        const Real longitudeOfAscendingNode
            = std::atan2( inclinationYComponent, inclinationXComponent );

        // Generate mean anomaly [rad].
        const Real meanAnomaly = meanAnomalyGenerator( randomSeed );

        // Compute true anomaly from mean anomaly [rad].
        const Real eccentricAnomaly = kep_toolbox::m2e( meanAnomaly, eccentricity );
        const Real trueAnomaly
            = astro::convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, eccentricity );

        query.bind( ":simulated",                       0 );
        query.bind( ":semi_major_axis",                 semiMajorAxis );
        query.bind( ":eccentricity",                    eccentricity );
        query.bind( ":inclination",                     inclination );
        query.bind( ":argument_of_periapsis",           argumentOfPeriapsis );
        query.bind( ":longitude_of_ascending_node",     longitudeOfAscendingNode );
        query.bind( ":true_anomaly",                    trueAnomaly );

        // Execute insert query.
        query.executeStep( );

        // Reset SQL insert query.
        query.reset( );
    }

        // Commit transaction.
        initialStatesTableInsertTransaction.commit( );
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

    const Real eccentricityFullWidthHalfMaximum
        = find( config, "eccentricity_full_width_half_maximum" )->value.GetDouble( );
    std::cout << "Eccentricity FWHM                  "
              << eccentricityFullWidthHalfMaximum << " [-]" << std::endl;

    const Real inclinationFullWidthHalfMaximum
        = find( config, "inclination_full_width_half_maximum" )->value.GetDouble( );
    std::cout << "Inclination FWHM                   "
              << inclinationFullWidthHalfMaximum << " [rad]" << std::endl;

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
                                       eccentricityFullWidthHalfMaximum,
                                       inclinationFullWidthHalfMaximum,
                                       integrator,
                                       startEpoch,
                                       endEpoch,
                                       timeStep,
                                       relativeTolerance,
                                       absoluteTolerance,
                                       databaseFilePath );
}

} // namespace dustsim