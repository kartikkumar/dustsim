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
     * @param[in] aNumberOfThreads            Number of parallel threads
     * @param[in] aNumberOfParticles          Number of dust particles to simulate
     * @param[in] aGravitationalParameter     Gravitational parameter of central body    [km^3 s^-2]
     * @param[in] aJ2AccelerationFlag         Flag indicating if J2 acceleration model is active
     * @param[in] aJ2Coefficient              J2 coefficient of gravity expansion                [-]
     * @param[in] anEquatorialRadius          Equatiorial radius for gravity expansion          [km]
     * @param[in] aRadiationPressureFlag      Flag indicating if radiation pressure acceleration
     *                                        model is active
     * @param[in] aParticleRadius             Radius of dust particle                       [micron]
     * @param[in] aParticleBulkDensity        Bulk density of dust particle                [kg m^-3]
     * @param[in] aRadiationPressureCoefficient
     *                                        Radiation pressure coefficient                     [-]
     * @param[in] aSemiMajorAxisMinimum       Minimum semi-major axis for uniform distribution  [km]
     * @param[in] aSemiMajorAxisMaximum       Maximum semi-major axis for uniform distribution  [km]
     * @param[in] anEccentricityFWHM          Full-Width Half-Maximum for normal distribution of
                                              eccentricity vector components                     [-]
     * @param[in] anInclinationFWHM           Full-Width Half-Maximum for normal distribution of
                                              inclination vector components                    [rad]
     * @param[in] anIntegrator                Name of selected numerical integrator
     * @param[in] aStartEpoch                 Start epoch for integration                        [s]
     * @param[in] aStepSize                   Fixed step size to generate integration output     [s]
     * @param[in] aNumberOfOutputSteps        Number of output steps requested
     * @param[in] aRelativeTolerance          Relative tolerance for integrator                  [-]
     * @param[in] anAbsoluteTolerance         Absolute tolerance for integrator                  [-]
     * @param[in] aDatabaseFilePath           Path to SQLite database for simulation results
     */
    BulkParticleSimulatorInput( const Real            aNumberOfThreads,
                                const Real            aNumberOfParticles,
                                const Real            aGravitationalParameter,
                                const bool            aJ2AccelerationFlag,
                                const Real            aJ2Coefficient,
                                const Real            anEquatorialRadius,
                                const bool            aRadiationPressureFlag,
                                const Real            aParticleRadius,
                                const Real            aParticleBulkDensity,
                                const Real            aRadiationPressureCoefficient,
                                const Real            aSemiMajorAxisMinimum,
                                const Real            aSemiMajorAxisMaximum,
                                const Real            anEccentricityFWHM,
                                const Real            anInclinationFWHM,
                                const Integrator      anIntegrator,
                                const Real            aStartEpoch,
                                const Real            aStepSize,
                                const Real            aNumberOfOutputSteps,
                                const Real            aRelativeTolerance,
                                const Real            anAbsoluteTolerance,
                                const std::string&    aDatabaseFilePath )
        : numberOfThreads( aNumberOfThreads ),
          numberOfParticles( aNumberOfParticles ),
          gravitationalParameter( aGravitationalParameter ),
          isJ2AccelerationModelActive( aJ2AccelerationFlag ),
          j2Coefficient( aJ2Coefficient ),
          equatorialRadius( anEquatorialRadius ),
          isRadiationPressureAccelerationModelActive( aRadiationPressureFlag ),
          particleRadius( aParticleRadius ),
          particleBulkDensity( aParticleBulkDensity ),
          radiationPressureCoefficient( aRadiationPressureCoefficient ),
          semiMajorAxisMinimum( aSemiMajorAxisMinimum ),
          semiMajorAxisMaximum( aSemiMajorAxisMaximum ),
          eccentricityFullWidthHalfMaximum( anEccentricityFWHM ),
          inclinationFullWidthHalfMaximum( anInclinationFWHM ),
          integrator( anIntegrator ),
          startEpoch( aStartEpoch ),
          stepSize( aStepSize ),
          outputSteps( aNumberOfOutputSteps ),
          relativeTolerance( aRelativeTolerance ),
          absoluteTolerance( anAbsoluteTolerance ),
          databaseFilePath( aDatabaseFilePath )
    { }

    //! Number of threads to parallelize simulations using OpenMP.
    const Int numberOfThreads;

    //! Number of dust particles to simulate.
    const Int numberOfParticles;

    //! Gravitational parameter of central body [km^3 s^-2].
    const Real gravitationalParameter;

    //! Boolean flag indicating if J2 acceleration model is active (true) or not (false).
    const bool isJ2AccelerationModelActive;

    //! J2-coefficient (unnormalized) of spherical harmonics expansion of gravity field [-].
    const Real j2Coefficient;

    //! Equatorial radius of central body corresponding with spherical harmonics gravity field [km].
    const Real equatorialRadius;

    //! Boolean flag indicating if radiation pressure acceleration model is active (true) or not
    //! (false).
    const bool isRadiationPressureAccelerationModelActive;

    //! Radius of dust particle [micron].
    const Real particleRadius;

    //! Bulk density of dust particle [kg m^-3].
    const Real particleBulkDensity;

    //! Radiation pressure coefficient [-].
    const Real radiationPressureCoefficient;

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

    //! Fixed step size to generate integration output [s].
    //! (N.B.: this is NOT the internal time step taken by the integration scheme but rather the
    //! step size between requested output points ).
    const Real stepSize;

    //! Number of equally-spaced output steps requested between start and end epoch.
    const Int outputSteps;

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
