/*
 * Copyright (c) 2017, K. Kumar, Delft University of Technology (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef DUSTSIM_BULK_PARTICLE_SIMULATOR_HPP
#define DUSTSIM_BULK_PARTICLE_SIMULATOR_HPP

#include <string>

#include <rapidjson/document.h>

#include "dustsim/typedefs.hpp"

namespace dustsim
{

//! Execute bulk particle simulator.
/*!
 * Executes simulations for a cloud of dust particle.
 *
 * This function is called when the user specifies the application mode to be
 * "bulk_particle_simulator".
 *
 * @param[in] config User-defined configuration options (extracted from JSON input file)
 */
void executeBulkParticleSimulator( const rapidjson::Document& config );

//! Input for bulk_particle_simulator application mode.
/*!
 * Data struct containing all valid input parameters to execute the bulk_particle_simulator
 * application mode. This struct is populated by the checkBulkParticleSimulatorInput() function.
 *
 * @sa checkBulkParticleSimulatorInput, executeBulkParticleSimulator
 */
struct BulkParticleSimulatorInput
{
public:

    //! Construct data struct.
    /*!
     * Constructs data struct based on verified input parameters.
     *
     * @sa checkSingleParticleSimulatorInput, executeSingleParticleSimulator
     * @param[in] aGravitationalParameter     Gravitational parameter of central body    [km^3 s^-2]
     * @param[in] aJ2AccelerationModelFlag    Flag indicating if J2 acceleration model is active
     * @param[in] aJ2Coefficient              J2 coefficient of gravity expansion                [-]
     * @param[in] anEquatorialRadius          Equatiorial radius for gravity expansion          [km]
     * @param[in] aNumberOfParticles          Number of dust particles to simulate
     * @param[in] aSemiMajorAxisMinimum       Minimum semi-major axis for uniform distribution  [km]
     * @param[in] aSemiMajorAxisMaximum       Maximum semi-major axis for uniform distribution  [km]
     * @param[in] anEccentricityFWHM          Full-Width Half-Maximum for normal distribution of
                                              eccentricity vector components                     [-]
     * @param[in] anInclinationFWHM           Full-Width Half-Maximum for normal distribution of
                                              inclination vector components                    [rad]
     * @param[in] anIntegrator                Name of selected numerical integrator
     * @param[in] aStartEpoch                 Start epoch for integration                        [s]
     * @param[in] anEndEpoch                  End epoch for integration                          [s]
     * @param[in] aTimeStep                   Time step for integration                          [s]
     * @param[in] aRelativeTolerance          Relative tolerance for integrator                  [-]
     * @param[in] anAbsoluteTolerance         Absolute tolerance for integrator                  [-]
     * @param[in] aDatabaseFilePath           Path to SQLite database for simulation results
     */
    BulkParticleSimulatorInput( const Real            aGravitationalParameter,
                                const bool            aJ2AccelerationModelFlag,
                                const Real            aJ2Coefficient,
                                const Real            anEquatorialRadius,
                                const Real            aNumberOfParticles,
                                const Real            aSemiMajorAxisMinimum,
                                const Real            aSemiMajorAxisMaximum,
                                const Real            anEccentricityFWHM,
                                const Real            anInclinationFWHM,
                                const Integrator      anIntegrator,
                                const Real            aStartEpoch,
                                const Real            anEndEpoch,
                                const Real            aTimeStep,
                                const Real            aRelativeTolerance,
                                const Real            anAbsoluteTolerance,
                                const std::string&    aDatabaseFilePath )
        : gravitationalParameter( aGravitationalParameter ),
          isJ2AccelerationModelActive( aJ2AccelerationModelFlag ),
          j2Coefficient( aJ2Coefficient ),
          equatorialRadius( anEquatorialRadius ),
          numberOfParticles( aNumberOfParticles ),
          semiMajorAxisMinimum( aSemiMajorAxisMinimum ),
          semiMajorAxisMaximum( aSemiMajorAxisMaximum ),
          eccentricityFullWidthHalfMaximum( anEccentricityFWHM ),
          inclinationFullWidthHalfMaximum( anInclinationFWHM ),
          integrator( anIntegrator ),
          startEpoch( aStartEpoch ),
          endEpoch( anEndEpoch ),
          timeStep( aTimeStep ),
          relativeTolerance( aRelativeTolerance ),
          absoluteTolerance( anAbsoluteTolerance ),
          databaseFilePath( aDatabaseFilePath )
    { }

    //! Gravitational parameter of central body [km^3 s^-2].
    const Real gravitationalParameter;

    //! Boolean flag indicating if J2 acceleration model is active (true) or not (false).
    const bool isJ2AccelerationModelActive;

    //! J2-coefficient (unnormalized) of spherical harmonics expansion of gravity field [-].
    const Real j2Coefficient;

    //! Equatorial radius of central body corresponding with spherical harmonics gravity field [km].
    const Real equatorialRadius;

    //! Number of dust particles to simulate.
    const Int numberOfParticles;

    //! Minimum semi-major axis corresponding to upper limit of uniform distribution [km].
    const Real semiMajorAxisMinimum;

    //! Maximum semi-major axis corresponding to upper limit of uniform distribution [km].
    const Real semiMajorAxisMaximum;

    //! Full-Width Half-Maximum for normal distribution of eccentricity vector components [-].
    const Real eccentricityFullWidthHalfMaximum;

    //! Full-Width Half-Maximum for normal distribution of inclination vector components [-].
    const Real inclinationFullWidthHalfMaximum;

    //! Selected numerical integrator.
    const Integrator integrator;

    //! Start epoch for simulator [s].
    const Real startEpoch;

    //! End epoch for simulator [s].
    const Real endEpoch;

    //! Time step for simulator [s].
    const Real timeStep;

    //! Relative tolerance for numerical integrator [-].
    const Real relativeTolerance;

    //! Absolute tolerance for numerical integrator [-].
    const Real absoluteTolerance;

    //! SQLite database file path.
    const std::string databaseFilePath;

protected:

private:
};

//! Check input parameters for bulk_particle_simulator application mode.
/*!
 * Checks that all inputs to execute bulk dust particle simulations are valid. If not, an error
 * is thrown with a short description of the problem. If all inputs are valid, a data struct
 * containing all the inputs is returned.
 *
 * @sa executeBulkParticleSimulator, BulkParticleSimulatorInput
 * @param[in] config User-defined configuration options (extracted from JSON input file)
 * @return           Struct containing all valid input for bulk_particle_simulator application
 *                   mode
 */
BulkParticleSimulatorInput checkBulkParticleSimulatorInput( const rapidjson::Document& config );

} // namespace dustsim

#endif // DUSTSIM_BULK_PARTICLE_SIMULATOR_HPP
