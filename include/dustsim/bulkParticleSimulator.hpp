/*
 * Copyright (c) 2009-2022 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#pragma once

#include <string>

#include <nlohmann/json.hpp>

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
 * @param[in]  config  User-defined configuration options (extracted from JSON input file)
 */
void executeBulkParticleSimulator(const nlohmann::json& config);

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
     * @param[in]  aNumberOfThreads             Number of parallel threads
     * @param[in]  aNumberOfParticles           Number of dust particles to simulate
     * @param[in]  aGravitationalParameter      Gravitational parameter of central body  [km^3 s^-2]
     * @param[in]  aJ2AccelerationFlag          Flag indicating if J2 acceleration model is active
     * @param[in]  aJ2Coefficient               J2 coefficient of gravity expansion              [-]
     * @param[in]  anEquatorialRadius           Equatiorial radius for gravity expansion         km]
     * @param[in]  aRadiationPressureFlag       Flag indicating if radiation pressure acceleration
     *                                          model is active
     * @param[in]  aParticleRadius              Radius of dust particle                     [micron]
     * @param[in]  aParticleBulkDensity         Bulk density of dust particle              [kg m^-3]
     * @param[in]  aRadiationPressureCoefficient
     *                                          Radiation pressure coefficient                   [-]
     * @param[in]  aSolarDistance               Average distance of central body from the Sun   [AU]
     * @param[in]  aSolarGravitationalParameter
     *                                          Gravitational parameter of the Sun       [km^3 s^-2]
     * @param[in]  aSolarEnergyFlux             Average energy flux at solar distance       [W m^-2]
     * @param[in]  aSemiMajorAxisMinimum        Minimum semi-major axis for uniform
     *                                          distribution                                    [km]
     * @param[in]  aSemiMajorAxisMaximum        Maximum semi-major axis for uniform
     *                                          distribution                                    [km]
     * @param[in]  anEccentricityFWHM           Full-Width Half-Maximum for normal distribution of
                                                eccentricity vector components                   [-]
     * @param[in]  anInclinationFWHM            Full-Width Half-Maximum for normal distribution of
                                                inclination vector components                  [rad]
     * @param[in]  anIntegrator                 Name of selected numerical integrator
     * @param[in]  aStartTime                   Start time for integration                       [s]
     * @param[in]  aTimeStep                    Time step for numerical integrator               [s]
     * @param[in]  anEndTime                    End time for integration                         [s]
     * @param[in]  aRelativeTolerance           Relative tolerance for integrator                [-]
     * @param[in]  anAbsoluteTolerance          Absolute tolerance for integrator                [-]
     * @param[in]  aMinimumStepSize             Minimum allowable step size for integrator       [s]
     * @param[in]  aMaximumStepSize             Maximum allowable step size for integrator       [s]
     * @param[in]  anOutputInterval             Interval at which output is written to database  [s]
     * @param[in]  aDatabaseFilePath            Path to SQLite database for simulation results
     */
    BulkParticleSimulatorInput(const Real            aNumberOfThreads,
                               const Real            aNumberOfParticles,
                               const Real            aGravitationalParameter,
                               const bool            aJ2AccelerationFlag,
                               const Real            aJ2Coefficient,
                               const Real            anEquatorialRadius,
                               const bool            aRadiationPressureFlag,
                               const Real            aParticleRadius,
                               const Real            aParticleBulkDensity,
                               const Real            aRadiationPressureCoefficient,
                               const Real            aSolarDistance,
                               const Real            aSolarGravitationalParameter,
                               const Real            aSolarEnergyFlux,
                               const Real            aSemiMajorAxisMinimum,
                               const Real            aSemiMajorAxisMaximum,
                               const Real            anEccentricityFWHM,
                               const Real            anInclinationFWHM,
                               const Integrator      anIntegrator,
                               const Real            aStartTime,
                               const Real            aTimeStep,
                               const Real            anEndTime,
                               const Real            aRelativeTolerance,
                               const Real            anAbsoluteTolerance,
                               const Real            aMinimumStepSize,
                               const Real            aMaximumStepSize,
                               const Real            anOutputInterval,
                               const std::string&    aDatabaseFilePath)
        : numberOfThreads(aNumberOfThreads ),
          numberOfParticles(aNumberOfParticles ),
          gravitationalParameter(aGravitationalParameter ),
          isJ2AccelerationModelActive(aJ2AccelerationFlag),
          j2Coefficient(aJ2Coefficient ),
          equatorialRadius(anEquatorialRadius ),
          isRadiationPressureAccelerationModelActive(aRadiationPressureFlag),
          particleRadius(aParticleRadius ),
          particleBulkDensity(aParticleBulkDensity ),
          radiationPressureCoefficient(aRadiationPressureCoefficient ),
          solarDistance(aSolarDistance),
          solarGravitationalParameter(aSolarGravitationalParameter),
          solarEnergyFlux(aSolarEnergyFlux),
          semiMajorAxisMinimum(aSemiMajorAxisMinimum ),
          semiMajorAxisMaximum(aSemiMajorAxisMaximum ),
          eccentricityFullWidthHalfMaximum(anEccentricityFWHM ),
          inclinationFullWidthHalfMaximum(anInclinationFWHM ),
          integrator(anIntegrator ),
          startTime(aStartTime ),
          timeStep(aTimeStep ),
          endTime(anEndTime ),
          relativeTolerance(aRelativeTolerance ),
          absoluteTolerance(anAbsoluteTolerance ),
          minimumStepSize(aMinimumStepSize ),
          maximumStepSize(aMaximumStepSize ),
          outputInterval(anOutputInterval ),
          databaseFilePath(aDatabaseFilePath )
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

    //! Mean solar distance [AU].
    const Real solarDistance;

    //! Solar gravitational parameter [km^3 s^-2].
    const Real solarGravitationalParameter;

    //! Mean solar energy flux [W m^-2].
    const Real solarEnergyFlux;

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

    //! Start time for simulator [s].
    const Real startTime;

    //! Step size to for integration scheme [s].
    //! For fixed step size integrators, the step size is constant throughout.
    //! For variable step size integrators, the step size is adapted based on the tolerances.
    const Real timeStep;

    //! End time for simulator [s].
    const Real endTime;

    //! Relative tolerance for numerical integrator [-].
    const Real relativeTolerance;

    //! Absolute tolerance for numerical integrator [-].
    const Real absoluteTolerance;

    //! Minimum allowable step size for numerical integrator [-].
    const Real minimumStepSize;

    //! Maximum allowable step size for numerical integrator [-].
    const Real maximumStepSize;

    //! Interval at which output is written to the database [s].
    const Real outputInterval;

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
 * @param[in]  config  User-defined configuration options (extracted from JSON input file)
 * @return             Struct containing all valid input for bulk_particle_simulator
 *                     application mode
 */
BulkParticleSimulatorInput checkBulkParticleSimulatorInput(const nlohmann::json& config);

} // namespace dustsim
