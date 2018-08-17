/*
 * Copyright (c) 2009-2018, K. Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <iomanip>
#include <iostream>
#include <iterator>
#include <cmath>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

#include <keplerian_toolbox.h>

#include <SQLiteCpp/SQLiteCpp.h>

#include <astro/astro.hpp>

#include <sml/sml.hpp>

#include <integrate/integrate.hpp>

#include "dustsim/bulkParticleSimulator.hpp"
#include "dustsim/dynamicalSystem.hpp"
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

    std::cout << "Setting up random number generators for initial states ..." << std::endl;

    // Initialize random seed.
    std::random_device randomNumberSeed;
    std::mt19937 randomNumberGenerator( randomNumberSeed( ) );

    // Set up generators for initial states for dust particles.
    // Semi-major axis values are generated from a uniform distribution, with minimum and maximum
    // set in the config file.
    std::uniform_real_distribution< Real > semiMajorAxisGenerator( input.semiMajorAxisMinimum,
                                                                   input.semiMajorAxisMaximum );

    // Eccentricity & argument of periapses values are generated from normal distributions for the
    // components of the eccentricity vector [h_e = e*cos( AoP ), k_e = e*sin( AoP )].
    const Real eccentricityStandardDeviation
        = input.eccentricityFullWidthHalfMaximum / ( 2.0 * std::sqrt( 2.0 * std::log( 2.0 ) ) );
    std::cout << "Eccentricity standard deviation    " << eccentricityStandardDeviation
              << " [-]" << std::endl;
    std::normal_distribution< Real > eccentricityComponentGenerator(
        0.0, eccentricityStandardDeviation );

    // Inclination & longitude of ascending node values are generated from normal distributions for
    // the components of the inclination vector [h_i = i*cos( RAAN ), k_i = i*sin( RAAN ) ].
    const Real inclinationStandardDeviation
        = input.inclinationFullWidthHalfMaximum / ( 2.0 * std::sqrt( 2.0 * std::log( 2.0 ) ) );
    std::cout << "Inclination standard deviation     " << inclinationStandardDeviation
              << " [rad]" << std::endl;
    std::normal_distribution< Real > inclinationComponentGenerator(
        0.0, inclinationStandardDeviation );

    // Mean anomaly values are generated from a uniform distribution.
    std::uniform_real_distribution< Real > meanAnomalyGenerator( 0.0, 2.0 * sml::SML_PI );

    std::cout << "Random number generators set up successfully!" << std::endl;
    std::cout << std::endl;

    std::cout << "Setting up database ..." << std::endl;

    // Open database in read/write mode.
    SQLite::Database database( input.databaseFilePath, SQLite::OPEN_READWRITE|SQLite::OPEN_CREATE );

    // Set names of tables in database.
    std::string metadataTableName = "metadata";
    std::string initialStatesTableName = "initial_states";
    std::string simulationResultsTableName = "simulation_results";

    std::cout << "Creating table '" << metadataTableName << "' ..." << std::endl;

    std::ostringstream metadataDropTable;
    metadataDropTable << "DROP TABLE IF EXISTS " << metadataTableName << ";";
    database.exec( metadataDropTable.str( ) );

    // Create metadata table.
    std::ostringstream metadataTableCreate;
    metadataTableCreate
        << "CREATE TABLE " << metadataTableName << " ("
        << "\"case_id\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
        << "\"gravitational_parameter\" REAL NOT NULL,"
        << "\"j2_coefficient\" REAL NOT NULL,"
        << "\"equatorial_radius\" REAL NOT NULL,"
        << "\"semi_major_axis_minimum\" REAL NOT NULL,"
        << "\"semi_major_axis_maximum\" REAL NOT NULL,"
        << "\"eccentricity_full_width_half_maximum\" REAL NOT NULL,"
        << "\"inclination_full_width_half_maximum\" REAL NOT NULL,"
        << "\"integrator\" TEXT NOT NULL,"
        << "\"start_time\" REAL NOT NULL,"
        << "\"initial_time_step\" REAL NOT NULL,"
        << "\"end_time\" INTEGER NOT NULL,"
        << "\"relative_tolerance\" REAL NOT NULL,"
        << "\"absolute_tolerance\" REAL NOT NULL,"
        << "\"minimum_step_size\" REAL NOT NULL,"
        << "\"maximum_step_size\" REAL NOT NULL);";

    database.exec( metadataTableCreate.str( ) );

    std::cout << "Creating table '" << initialStatesTableName << "' ..." << std::endl;

    std::ostringstream initialStatesDropTable;
    initialStatesDropTable << "DROP TABLE IF EXISTS " << initialStatesTableName << ";";
    database.exec( initialStatesDropTable.str( ) );

    // Create initial states table.
    std::ostringstream initialStatesTableCreate;
    initialStatesTableCreate
        << "CREATE TABLE " << initialStatesTableName << " ("
        << "\"simulation_id\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
        << "\"is_simulated\" INTEGER NOT NULL,"
        << "\"semi_major_axis\" REAL NOT NULL,"
        << "\"eccentricity\" REAL NOT NULL,"
        << "\"inclination\" REAL NOT NULL,"
        << "\"argument_of_periapsis\" REAL NOT NULL,"
        << "\"longitude_of_ascending_node\" REAL NOT NULL,"
        << "\"true_anomaly\" REAL NOT NULL);";

    database.exec( initialStatesTableCreate.str( ) );

    std::cout << "Creating table '" << simulationResultsTableName << "' ..." << std::endl;

    std::ostringstream simulationResultsDropTable;
    simulationResultsDropTable << "DROP TABLE IF EXISTS " << simulationResultsTableName << ";";
    database.exec( simulationResultsDropTable.str( ) );

    // Create simulation results table.
    std::ostringstream simulationResultsTableCreate;
    simulationResultsTableCreate
        << "CREATE TABLE " << simulationResultsTableName << " ("
        << "\"time_id\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
        << "\"simulation_id\" INTEGER NOT NULL,"
        << "\"time\" REAL NOT NULL,"
        << "\"semi_major_axis\" REAL NOT NULL,"
        << "\"eccentricity\" REAL NOT NULL,"
        << "\"inclination\" REAL NOT NULL,"
        << "\"argument_of_periapsis\" REAL NOT NULL,"
        << "\"longitude_of_ascending_node\" REAL NOT NULL,"
        << "\"true_anomaly\" REAL NOT NULL,"
        << "\"x_position\" REAL NOT NULL,"
        << "\"y_position\" REAL NOT NULL,"
        << "\"z_position\" REAL NOT NULL,"
        << "\"x_velocity\" REAL NOT NULL,"
        << "\"y_velocity\" REAL NOT NULL,"
        << "\"z_velocity\" REAL NOT NULL,"
        << "\"x_unit_position_sun\" REAL NOT NULL,"
        << "\"y_unit_position_sun\" REAL NOT NULL,"
        << "\"z_unit_position_sun\" REAL NOT NULL);";

    database.exec( simulationResultsTableCreate.str( ) );

    std::cout << "Database set up successfully!" << std::endl;
    std::cout << std::endl;

    std::cout << "Populating metadata table ..." << std::endl;

    std::ostringstream metadataInsertString;
    metadataInsertString
        << "INSERT INTO '" << metadataTableName << "' VALUES ("
        << "NULL,";
    metadataInsertString
        << std::setprecision( std::numeric_limits< double >::digits10 )
        << input.gravitationalParameter << ","
        << input.j2Coefficient << ","
        << input.equatorialRadius << ","
        << input.semiMajorAxisMinimum << ","
        << input.semiMajorAxisMaximum << ","
        << input.eccentricityFullWidthHalfMaximum << ","
        << input.inclinationFullWidthHalfMaximum << ",";
    metadataInsertString
        << input.integrator << ",";
    metadataInsertString
        << std::setprecision( std::numeric_limits< double >::digits10 )
        << input.startTime << ","
        << input.timeStep << ","
        << input.endTime << ","
        << input.relativeTolerance << ","
        << input.absoluteTolerance << ","
        << input.minimumStepSize << ","
        << input.maximumStepSize
        << ");";

    database.exec( metadataInsertString.str( ) );

    std::cout << "Populating initial states table ..." << std::endl;

    std::ostringstream initialStatesInsertString;
    initialStatesInsertString
        << "INSERT INTO '" << initialStatesTableName << "' VALUES ("
        << ":simulation_id,"
        << ":is_simulated,"
        << ":semi_major_axis,"
        << ":eccentricity,"
        << ":inclination,"
        << ":argument_of_periapsis,"
        << ":longitude_of_ascending_node,"
        << ":true_anomaly"
        << ");";

    SQLite::Statement initialStatesInsertQuery( database, initialStatesInsertString.str( ) );

    // Set up database transaction.
    SQLite::Transaction initialStatesInsertTransaction( database );

    // Create storage container to store initial states (Keplerian elements) generated.
    typedef std::map< Int, State > IntialStates;
    IntialStates initialStatesInKeplerianElements;

    for ( int i = 0; i < input.numberOfParticles ; ++i )
    {
        // Store simulation id.
        const Real simulationId = i+1;

        // Generate random semi-major axis [km].
        const Real semiMajorAxis = semiMajorAxisGenerator( randomNumberGenerator );

        // Generate random eccentricity vector components.
        const Real eccentricityXComponent = eccentricityComponentGenerator( randomNumberGenerator );
        const Real eccentricityYComponent = eccentricityComponentGenerator( randomNumberGenerator );

        // Compute eccentricity [-] and argument of periapsis [rad].
        const Real eccentricity = std::sqrt( eccentricityXComponent * eccentricityXComponent
                                             + eccentricityYComponent * eccentricityYComponent );

        const Real argumentOfPeriapsis
            = std::atan2( eccentricityYComponent, eccentricityXComponent );

        // Generate random inclination vector components.
        const Real inclinationXComponent = inclinationComponentGenerator( randomNumberGenerator );
        const Real inclinationYComponent = inclinationComponentGenerator( randomNumberGenerator );

        // Compute inclination [rad] and longitude of ascending node [rad].
        const Real inclination = std::sqrt( inclinationXComponent * inclinationXComponent
                                             + inclinationYComponent * inclinationYComponent );

        const Real longitudeOfAscendingNode
            = std::atan2( inclinationYComponent, inclinationXComponent );

        // Generate mean anomaly [rad].
        const Real meanAnomaly = meanAnomalyGenerator( randomNumberGenerator );

        // Compute true anomaly from mean anomaly [rad].
        // @TODO: replace mean-to-eccentric anomaly conversion with implementation in
        // openastro/astro.
        const Real eccentricAnomaly = kep_toolbox::m2e( meanAnomaly, eccentricity );
        const Real trueAnomaly
            = astro::convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, eccentricity );

        // Store initial state in container.
        initialStatesInKeplerianElements.insert( { simulationId, State( { semiMajorAxis,
                                                                          eccentricity,
                                                                          inclination,
                                                                          argumentOfPeriapsis,
                                                                          longitudeOfAscendingNode,
                                                                          trueAnomaly } )  } );

        // Bind initial state values to database insert query.
        initialStatesInsertQuery.bind( ":simulation_id",               simulationId );
        initialStatesInsertQuery.bind( ":is_simulated",                0 );
        initialStatesInsertQuery.bind( ":semi_major_axis",             semiMajorAxis );
        initialStatesInsertQuery.bind( ":eccentricity",                eccentricity );
        initialStatesInsertQuery.bind( ":inclination",                 inclination );
        initialStatesInsertQuery.bind( ":argument_of_periapsis",       argumentOfPeriapsis );
        initialStatesInsertQuery.bind( ":longitude_of_ascending_node", longitudeOfAscendingNode );
        initialStatesInsertQuery.bind( ":true_anomaly",                trueAnomaly );

        // Execute insert query.
        initialStatesInsertQuery.executeStep( );

        // Reset SQL insert query.
        initialStatesInsertQuery.reset( );
    }

    // Commit transaction.
    initialStatesInsertTransaction.commit( );

    std::cout << "Initial states table populated successfully!" << std::endl;
    std::cout << std::endl;

    std::cout << "Running simulations ..." << std::endl;

    std::ostringstream simulationResultsInsertString;
    simulationResultsInsertString
        << "INSERT INTO '" << simulationResultsTableName << "' VALUES ("
        << "NULL,"
        << ":simulation_id,"
        << ":time,"
        << ":semi_major_axis,"
        << ":eccentricity,"
        << ":inclination,"
        << ":argument_of_periapsis,"
        << ":longitude_of_ascending_node,"
        << ":true_anomaly,"
        << ":x_position,"
        << ":y_position,"
        << ":z_position,"
        << ":x_velocity,"
        << ":y_velocity,"
        << ":z_velocity,"
        << ":x_unit_position_sun,"
        << ":y_unit_position_sun,"
        << ":z_unit_position_sun"
        << ");";

    SQLite::Statement simulationResultsInsertQuery(
        database, simulationResultsInsertString.str( ) );

    // Set up database transaction.
    SQLite::Transaction simulationResultsInsertTransaction( database );

    // Set up query to update completed states in initial states table once simulation has been run.
    std::ostringstream initialStatesUpdateString;
    initialStatesUpdateString << "UPDATE " << initialStatesTableName
                              << " SET \"is_simulated\" = 1 WHERE simulation_id = :simulation_id;";

    // Set up database query.
    SQLite::Statement initialStatesUpdateQuery( database, initialStatesUpdateString.str( ) );

#pragma omp parallel for schedule( dynamic ) num_threads( input.numberOfThreads )
    // Loop through the table retrieved from the database, step-by-step and execute simulations.
    // OpenMP support for range-based for loops is not well-established, so using integer loop
    // instead.
    for ( unsigned int j = 0; j < initialStatesInKeplerianElements.size( ); ++j )
    {
        IntialStates::iterator iterator = initialStatesInKeplerianElements.begin( );
        std::advance( iterator, j );

        const Int simulationId = iterator->first;
        const State initialStateInKeplerianElements = iterator->second;

#pragma omp critical( outputToConsole )
        {
            std::cout << "Executing simulation ID: " << simulationId << std::endl;
        }

        const State initialState = astro::convertKeplerianToCartesianElements(
            initialStateInKeplerianElements, input.gravitationalParameter );

        // Set current time and state and step size to initial values specified by user.
        State state = initialState;
        Real time = input.startTime;
        Real stepSize = input.timeStep;

        auto stateDerivativePointer = std::bind( &DynamicalSystem::operator( ),
                                                 &dynamics,
                                                 std::placeholders::_1,
                                                 std::placeholders::_2 );

        Int outputIntervalCounter = 1;

        while ( time < input.endTime )
        {
            Real previousTime = time;
            State previousState = state;

            integrate::stepRKF78< Real, State >( time,
                                                 state,
                                                 stepSize,
                                                 stateDerivativePointer,
                                                 input.relativeTolerance,
                                                 input.minimumStepSize,
                                                 input.maximumStepSize );

            const Real nextOutputTime = input.outputInterval * outputIntervalCounter;

            if ( time > nextOutputTime + input.minimumStepSize )
            {
                Real outputStepSize = ( nextOutputTime - previousTime );
                integrate::stepRKF78< Real, State >( previousTime,
                                                     previousState,
                                                     outputStepSize,
                                                     stateDerivativePointer,
                                                     input.relativeTolerance,
                                                     input.minimumStepSize,
                                                     input.maximumStepSize );

                const State stateInKeplerianElements
                    = astro::convertCartesianToKeplerianElements(
                        previousState, input.gravitationalParameter );

                // To avoid locking of the database, this section is thread-critical, so will be
                // executed one-by-one by multiple threads.
#pragma omp critical( writeOutputToDatabase )
                {
                    simulationResultsInsertQuery.bind( ":simulation_id", simulationId );
                    simulationResultsInsertQuery.bind( ":time", previousTime );
                    simulationResultsInsertQuery.bind(
                        ":semi_major_axis",
                        stateInKeplerianElements[ astro::semiMajorAxisIndex ] );
                    simulationResultsInsertQuery.bind(
                        ":eccentricity",
                        stateInKeplerianElements[ astro::eccentricityIndex ] );
                    simulationResultsInsertQuery.bind(
                        ":inclination",
                        stateInKeplerianElements[ astro::inclinationIndex ] );
                    simulationResultsInsertQuery.bind(
                        ":argument_of_periapsis",
                        stateInKeplerianElements[ astro::argumentOfPeriapsisIndex ] );
                    simulationResultsInsertQuery.bind(
                        ":longitude_of_ascending_node",
                        stateInKeplerianElements[ astro::longitudeOfAscendingNodeIndex ] );
                    simulationResultsInsertQuery.bind(
                        ":true_anomaly",
                        stateInKeplerianElements[ astro::trueAnomalyIndex ] );
                    simulationResultsInsertQuery.bind(
                        ":x_position", state[ astro::xPositionIndex ] );
                    simulationResultsInsertQuery.bind(
                        ":y_position", state[ astro::yPositionIndex ] );
                    simulationResultsInsertQuery.bind(
                        ":z_position", state[ astro::zPositionIndex ] );
                    simulationResultsInsertQuery.bind(
                        ":x_velocity", state[ astro::xVelocityIndex ] );
                    simulationResultsInsertQuery.bind(
                        ":y_velocity", state[ astro::yVelocityIndex ] );
                    simulationResultsInsertQuery.bind(
                        ":z_velocity", state[ astro::zVelocityIndex ] );
                    simulationResultsInsertQuery.bind( ":x_unit_position_sun", 0.0 );
                    simulationResultsInsertQuery.bind( ":y_unit_position_sun", 0.0 );
                    simulationResultsInsertQuery.bind( ":z_unit_position_sun", 0.0 );

                    // Execute insert query.
                    simulationResultsInsertQuery.executeStep( );

                    // Reset SQL insert query.
                    simulationResultsInsertQuery.reset( );

                    initialStatesUpdateQuery.bind( ":simulation_id", simulationId );

                    // Execute insert query.
                    initialStatesUpdateQuery.executeStep( );

                    // Reset SQL insert query.
                    initialStatesUpdateQuery.reset( );
                }

                outputIntervalCounter++;
            }
        }

    }

    // Commit transaction.
    simulationResultsInsertTransaction.commit( );

    std::cout << "Simulations run successfully!" << std::endl;
}

//! Check input parameters for bulk_particle_simulator application mode.
BulkParticleSimulatorInput checkBulkParticleSimulatorInput( const rapidjson::Document& config )
{
    // Extract number of threads to parallelize simulations using OpenMP.
    const Int numberOfThreads = find( config, "threads" )->value.GetInt( );
    std::cout << "Number of threads                  " << numberOfThreads << std::endl;

    // Extract number of particles to simulate.
    const Int numberOfParticles = find( config, "particles" )->value.GetInt( );
    std::cout << "Number of particles                " << numberOfParticles << std::endl;

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
        if ( integratorString.compare( "rkf45" ) == 0 )
        {
            integrator = rkf45;
        }
        else if ( integratorString.compare( "rkf78" ) == 0 )
        {
            integrator = rkf78;
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
    const Real startTime            = find( config, "start_time" )->value.GetDouble( );
    std::cout << "Start time                         " << startTime << " [s]" << std::endl;
    const Real endTime              = find( config, "end_time" )->value.GetDouble( );
    std::cout << "End time                           " << endTime << " [s]" << std::endl;
    const Real timeStep             = find( config, "time_step" )->value.GetDouble( );
    std::cout << "Time step                          " << timeStep << " [s]" << std::endl;

    // Extract integrator tolerances.
    const Real relativeTolerance    = find( config, "relative_tolerance" )->value.GetDouble( );
    std::cout << "Relative tolerance                 " << relativeTolerance << " [-]" << std::endl;
    const Real absoluteTolerance    = find( config, "absolute_tolerance" )->value.GetDouble( );
    std::cout << "Absolute tolerance                 " << absoluteTolerance << " [-]" << std::endl;

    // Extract integrator step size bounds.
    const Real minimumStepSize    = find( config, "minimum_step_size" )->value.GetDouble( );
    std::cout << "Minimum step size                  " << minimumStepSize << " [-]" << std::endl;
    const Real maximumStepSize    = find( config, "maximum_step_size" )->value.GetDouble( );
    std::cout << "Maximum step size                  " << maximumStepSize << " [-]" << std::endl;

    const Real outputInterval     = find( config, "output_interval" )->value.GetDouble( );
    std::cout << "Output interval                    " << outputInterval  << " [s]" << std::endl;

    // Extract SQLite database settings.
    const std::string databaseFilePath = find( config, "database_file_path" )->value.GetString( );
    std::cout << "Database file path                 " << databaseFilePath << std::endl;

    return BulkParticleSimulatorInput( numberOfThreads,
                                       numberOfParticles,
                                       gravitationalParameter,
                                       j2AcclerationModelFlag,
                                       j2Coefficient,
                                       equatorialRadius,
                                       radiationPressureFlag,
                                       particleRadius,
                                       particleBulkDensity,
                                       radiationPressureCoefficient,
                                       semiMajorAxisMinimum,
                                       semiMajorAxisMaximum,
                                       eccentricityFullWidthHalfMaximum,
                                       inclinationFullWidthHalfMaximum,
                                       integrator,
                                       startTime,
                                       timeStep,
                                       endTime,
                                       relativeTolerance,
                                       absoluteTolerance,
                                       minimumStepSize,
                                       maximumStepSize,
                                       outputInterval,
                                       databaseFilePath );
}

} // namespace dustsim
