/*
 * Copyright (c) 2017, K. Kumar, Delft University of Technology (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <keplerian_toolbox.h>

#include <SQLiteCpp/SQLiteCpp.h>

#include <astro/astro.hpp>

#include <sml/sml.hpp>

#include "dustsim/bulkParticleSimulator.hpp"
#include "dustsim/dynamicalSystem.hpp"
#include "dustsim/tools.hpp"
#include "dustsim/outputWriter.hpp"

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
        << "\"start_epoch\" REAL NOT NULL,"
        << "\"step_size\" REAL NOT NULL,"
        << "\"output_steps\" REAL NOT NULL,"
        << "\"relative_tolerance\" REAL NOT NULL,"
        << "\"absolute_tolerance\" REAL NOT NULL);";

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
        << "\"simulated\" INTEGER NOT NULL,"
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
        << "\"ephemeris_id\" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,"
        << "\"simulation_id\" INTEGER NOT NULL,"
        << "\"epoch\" REAL NOT NULL,"
        << "\"semi_major_axis\" REAL NOT NULL,"
        << "\"eccentricity\" REAL NOT NULL,"
        << "\"inclination\" REAL NOT NULL,"
        << "\"argument_of_periapsis\" REAL NOT NULL,"
        << "\"longitude_of_ascending_node\" REAL NOT NULL,"
        << "\"true_anomaly\" REAL NOT NULL);";

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
        << input.startEpoch << ","
        << input.stepSize << ","
        << input.outputSteps << ","
        << input.relativeTolerance << ","
        << input.absoluteTolerance
        << ");";

    database.exec( metadataInsertString.str( ) );

    std::cout << "Populating initial states table ..." << std::endl;

    std::ostringstream initialStatesInsertString;
    initialStatesInsertString
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

    SQLite::Statement initialStatesInsertQuery( database, initialStatesInsertString.str( ) );

    // Set up database transaction.
    SQLite::Transaction initialStatesInsertTransaction( database );

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

        initialStatesInsertQuery.bind( ":simulated",                   0 );
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

    // Create instance of dynamical system.
    std::cout << "Set up dynamical model ..." << std::endl;
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

    std::cout << "Running simulations ..." << std::endl;

    std::ostringstream simulationResultsInsertString;
    simulationResultsInsertString
        << "INSERT INTO '" << simulationResultsTableName << "' VALUES ("
        << "NULL,"
        << ":simulation_id,"
        << ":epoch,"
        << ":semi_major_axis,"
        << ":eccentricity,"
        << ":inclination,"
        << ":argument_of_periapsis,"
        << ":longitude_of_ascending_node,"
        << ":true_anomaly"
        << ");";

    SQLite::Statement simulationResultsInsertQuery(
        database, simulationResultsInsertString.str( ) );

    // Set up database transaction.
    SQLite::Transaction simulationResultsInsertTransaction( database );

    // Set up query to fetch initial states table from database.
    std::ostringstream initialStatesFetchString;
    initialStatesFetchString << "SELECT * FROM " << initialStatesTableName << ";";

    // Set up database query.
    SQLite::Statement initialStatesFetchQuery( database, initialStatesFetchString.str( ) );

    // Fetch initial states table in database.
    typedef std::map< Int, State > InitialStates;
    InitialStates initialStates;
    while ( initialStatesFetchQuery.executeStep( ) )
    {
        const Int simulationId = initialStatesFetchQuery.getColumn( 0 );
        State initialStateKeplerianElements;
        initialStateKeplerianElements[ astro::semiMajorAxisIndex ]
            = initialStatesFetchQuery.getColumn( 2 );
        initialStateKeplerianElements[ astro::eccentricityIndex ]
            = initialStatesFetchQuery.getColumn( 3 );
        initialStateKeplerianElements[ astro::inclinationIndex ]
            = initialStatesFetchQuery.getColumn( 4 );
        initialStateKeplerianElements[ astro::argumentOfPeriapsisIndex ]
            = initialStatesFetchQuery.getColumn( 5 );
        initialStateKeplerianElements[ astro::longitudeOfAscendingNodeIndex ]
            = initialStatesFetchQuery.getColumn( 6 );
        initialStateKeplerianElements[ astro::trueAnomalyIndex ]
            = initialStatesFetchQuery.getColumn( 7 );

        initialStates[ simulationId ] = initialStateKeplerianElements;
    }

    InitialStates::iterator initialState = initialStates.begin( );

    #pragma omp parallel for num_threads( input.numberOfThreads )
    // Loop through the table retrieved from the database, step-by-step and execute simulations.
    for ( unsigned int j = 0; j < initialStates.size( ); ++j )
    {
        // Compute initial state in Cartesian elements.
        State currentState = astro::convertKeplerianToCartesianElements(
            initialState->second, input.gravitationalParameter );

        std::ostringstream integrationOutput;
        StateHistoryWriter writer( integrationOutput, input.gravitationalParameter );

        // Execute selected numerical integrator.
        #pragma omp critical( outputToConsole )
        {
            std::cout << "Executing numerical integration " << j
                      << ": ID " << initialState->first << std::endl;
        }

        if ( input.integrator == rk4 )
        {
            using namespace boost::numeric::odeint;
            integrate_n_steps( runge_kutta4< State >( ),
                               dynamics,
                               currentState,
                               input.startEpoch,
                               input.stepSize,
                               input.outputSteps,
                               writer );
        }
        else if ( input.integrator == dopri5 )
        {
            using namespace boost::numeric::odeint;
            integrate_n_steps( make_dense_output( input.relativeTolerance,
                                                  input.absoluteTolerance,
                                                  runge_kutta_dopri5< State >( ) ),
                               dynamics,
                               currentState,
                               input.startEpoch,
                               input.stepSize,
                               input.outputSteps,
                               writer );
        }
        else if ( input.integrator == bs )
        {
            using namespace boost::numeric::odeint;
            bulirsch_stoer_dense_out< State > stepper( input.absoluteTolerance,
                                                       input.relativeTolerance );
            integrate_n_steps( stepper,
                               dynamics,
                               currentState,
                               input.startEpoch,
                               input.stepSize,
                               input.outputSteps,
                               writer );
        }
        else
        {
            throw std::runtime_error( "Selected numerical integrator is invalid!" );
        }

        // Write output generated by numerical integrator to the simulation results table in the
        // database.
        const std::string output = integrationOutput.str( );
        std::vector< std::string > lines;
        boost::algorithm::split( lines, output, boost::is_any_of( "\n" ) );

        for ( unsigned int i = 0; i < lines.size( ) - 1; ++i )
        {
            std::vector< std::string > stateElements;
            boost::algorithm::split( stateElements, lines.at( i ), boost::is_any_of( "," ) );

            // To avoid locking of the database, this section is thread-critical, so will be
            // executed one-by-one by multiple threads.
            #pragma omp critical( writeOutputToDatabase )
            {
                simulationResultsInsertQuery.bind(
                        ":simulation_id",               initialState->first );
                simulationResultsInsertQuery.bind(
                        ":epoch",                       i * input.stepSize );
                simulationResultsInsertQuery.bind(
                        ":semi_major_axis",             stateElements.at( 7 ) );
                simulationResultsInsertQuery.bind(
                        ":eccentricity",                stateElements.at( 8 ) );
                simulationResultsInsertQuery.bind(
                        ":inclination",                 stateElements.at( 9 ) );
                simulationResultsInsertQuery.bind(
                        ":argument_of_periapsis",       stateElements.at( 10 ) );
                simulationResultsInsertQuery.bind(
                        ":longitude_of_ascending_node", stateElements.at( 11 ) );
                simulationResultsInsertQuery.bind(
                        ":true_anomaly",                stateElements.at( 12 ) );

                // Execute insert query.
                simulationResultsInsertQuery.executeStep( );

                // Reset SQL insert query.
                simulationResultsInsertQuery.reset( );
            }
        }

        initialState++;
    }

    // // Commit transaction.
    // simulationResultsInsertTransaction.commit( );

    // std::cout << "Simulations run successfully!" << std::endl;
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
    const Real stepSize             = find( config, "step_size" )->value.GetDouble( );
    std::cout << "Step size                          " << stepSize << " [s]" << std::endl;
    const Real outputSteps          = find( config, "output_steps" )->value.GetInt( );
    std::cout << "Number of output steps             " << outputSteps << std::endl;

    // Extract integrator tolerances.
    const Real relativeTolerance    = find( config, "relative_tolerance" )->value.GetDouble( );
    std::cout << "Relative tolerance                 " << relativeTolerance << " [-]" << std::endl;
    const Real absoluteTolerance    = find( config, "absolute_tolerance" )->value.GetDouble( );
    std::cout << "Absolute tolerance                 " << absoluteTolerance << " [-]" << std::endl;

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
                                       startEpoch,
                                       stepSize,
                                       outputSteps,
                                       relativeTolerance,
                                       absoluteTolerance,
                                       databaseFilePath );
}

} // namespace dustsim