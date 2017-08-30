/*
 * Copyright (c) 2017, K. Kumar, Delft University of Technology (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef DUSTSIM_SINGLE_PARTICLE_SIMULATOR_HPP
#define DUSTSIM_SINGLE_PARTICLE_SIMULATOR_HPP

#include <string>

#include <rapidjson/document.h>

#include "dustsim/typedefs.hpp"

namespace dustsim
{

//! Execute single particle simulator.
/*!
 * Executes a single dust particle simulation.
 *
 * This function is called when the user specifies the application mode to be
 * "single_particle_simulator".
 *
 * @param[in] config User-defined configuration options (extracted from JSON input file)
 */
void executeSingleParticleSimulator( const rapidjson::Document& config );

// //! Input for single_particle_simulator application mode.
// /*!
//  * Data struct containing all valid input parameters to execute the single_particle_simulator
//  * application mode. This struct is populated by the checkSingleParticleSimulatorInput() function.
//  *
//  * @sa checkSingleParticleSimulatorInput, executeSingleParticleSimulator
//  */

struct SingleParticleSimulatorInput
{
public:

    //! Construct data struct.
    /*!
     * Constructs data struct based on verified input parameters.
     *
     * @sa checkSingleParticleSimulatorInput, executeSingleParticleSimulator
     * @param[in] aGravitationalParameter     Gravitational parameter of central body    [km^3 s^-2]
     * @param[in] aJ2Coefficient              J2 coefficient of gravity expansion                [-]
     * @param[in] anEquatorialRadius          Equatiorial radius for gravity expansion          [km]
     * @param[in] anInitialKeplerState        Initial state in Keplerian elements
     * @param[in] aStartEpoch                 Start epoch for integration                        [s]
     * @param[in] anEndEpoch                  End epoch for integration                          [s]
     * @param[in] aTimeStep                   Time step for integration                          [s]
     * @param[in] aRelativeTolerance          Relative tolerance for integrator                  [-]
     * @param[in] anAbsoluteTolerance         Absolute tolerance for integrator                  [-]
     * @param[in] aStateHistoryFilePath       Path to output file for state history
     */
    SingleParticleSimulatorInput( const Real            aGravitationalParameter,
                                  const Real            aJ2Coefficient,
                                  const Real            anEquatorialRadius,
                                  const State&          anInitialKeplerState,
                                  const Real            aStartEpoch,
                                  const Real            anEndEpoch,
                                  const Real            aTimeStep,
                                  const Real            aRelativeTolerance,
                                  const Real            anAbsoluteTolerance,
                                  const std::string&    aStateHistoryFilePath )
        : gravitationalParameter( aGravitationalParameter ),
          j2Coefficient( aJ2Coefficient ),
          equatorialRadius( anEquatorialRadius ),
          initialStateKeplerianElements( anInitialKeplerState ),
          startEpoch( aStartEpoch ),
          endEpoch( anEndEpoch ),
          timeStep( aTimeStep ),
          relativeTolerance( aRelativeTolerance ),
          absoluteTolerance( anAbsoluteTolerance ),
          stateHistoryFilePath( aStateHistoryFilePath )
    { }

    //! Gravitational parameter of central body [km^3 s^-2].
    const Real gravitationalParameter;

    //! J2-coefficient (unnormalized) of spherical harmonics expansion of gravity field [-].
    const Real j2Coefficient;

    //! Equatorial radius of central body corresponding with spherical harmonics gravity field [km].
    const Real equatorialRadius;

    //! Initial state in Keplerian elements [km, -, rad, rad, rad, rad].
    const State initialStateKeplerianElements;

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

    //! State history file path.
    const std::string stateHistoryFilePath;

protected:

private:
};

//! Check input parameters for single_particle_simulator application mode.
/*!
 * Checks that all inputs to execute a single dust particle simulation are valid. If not, an error
 * is thrown with a short description of the problem. If all inputs are valid, a data struct
 * containing all the inputs is returned.
 *
 * @sa executeSingleParticleSimulator, SingleParticleSimulatorInput
 * @param[in] config User-defined configuration options (extracted from JSON input file)
 * @return           Struct containing all valid input for single_particle_simulator application
 *                   mode
 */
SingleParticleSimulatorInput checkSingleParticleSimulatorInput( const rapidjson::Document& config );

} // namespace dustsim

#endif // DUSTSIM_SINGLE_PARTICLE_SIMULATOR_HPP
