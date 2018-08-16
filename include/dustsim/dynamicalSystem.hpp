/*
 * Copyright (c) 2017, K. Kumar, Delft University of Technology (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef DUSTSIM_DYNAMICAL_SYSTEM_HPP
#define DUSTSIM_DYNAMICAL_SYSTEM_HPP

#include <astro/astro.hpp>

#include "dustsim/state.hpp"
#include "dustsim/typedefs.hpp"

namespace dustsim
{

//! Class containing parameters and models to describe the dynamical system.
/*!
 * This class contains parameters and models to describe the dynamical system that governs the
 * motion of dust particles. The dynamical system is defined in terms of a set of 1st-order ODEs,
 * with the right-hand side including all the force models acting on the dust particles. This
 * dynamical system is provided as input to a numerical integrator to compute the derivative of the
 * state vector.
 */
class DynamicalSystem
{
public:

    //! Construct dynamical system.
    /*!
     * Constructor for dynamical system, taking model parameters to define the dynamical
     * environment.
     *
     * @param[in] aGravitationalParameter   Gravitational parameter of the central body  [km^3 s^-2]
     * @param[in] aJ2AccelerationFlag       Flag indicating if J2 acceleration model is active
     * @param[in] aJ2Coefficient            J2 coefficient of gravity expansion                  [-]
     * @param[in] anEquatorialRadius        Equatiorial radius for gravity expansion            [km]
     * @param[in] aRadiationPressureFlag    Flag indicating if radiation pressure acceleration
     *                                      model is active
     * @param[in] aParticleRadius           Radius of dust particle                         [micron]
     * @param[in] aParticleBulkDensity      Bulk density of dust particle                   [kg m^-3]
     * @param[in] aRadiationPressureCoefficient
     *                                      Radiation pressure coefficient                       [-]
     */
    DynamicalSystem( const Real aGravitationalParameter,
                     const bool aJ2AccelerationFlag,
                     const Real aJ2Coefficient,
                     const Real anEquatorialRadius,
                     const bool aRadiationPressureFlag,
                     const Real aParticleRadius,
                     const Real aParticleBulkDensity,
                     const Real aRadiationPressureCoefficient )
        : gravitationalParameter( aGravitationalParameter ),
          isJ2AccelerationModelActive( aJ2AccelerationFlag ),
          j2Coefficient( aJ2Coefficient ),
          equatorialRadius( anEquatorialRadius ),
          isRadiationPressureAccelerationModelActive( aRadiationPressureFlag ),
          particleRadius( aParticleRadius ),
          particleBulkDensity( aParticleBulkDensity ),
          radiationPressureCoefficient( aRadiationPressureCoefficient )
    { }

    //! Overload ()-operator to compute state derivative using dynamical system.
    /*!
     * Overloads the ()-operator to compute the state derivative, given the current state and epoch,
     * based on the parameters and models used to construct the dynamical system.
     * This function fulfills the prototype for the numerical integrators in the openastro integrate
     * library.
     *
     * @param[in]  time              Current simulation epoch
     * @param[in]  state             Current state of the dynamical system (1-D vector)
     * @return                       Computed state derivative of the dynamical system (1-D vector)
     */
    State operator( )( const Real time, const State& state )
    {
        const State position( { state[ astro::xPositionIndex ],
                                state[ astro::yPositionIndex ],
                                state[ astro::zPositionIndex ] } );

        // Compute the total acceleration acting on the system as a sum of the forces.
        // Central body gravity is included by default.
        State acceleration( astro::computeCentralBodyAcceleration( gravitationalParameter,
                                                                   position ) );

        // Add J2 acceleration if model is set to active.
        if ( isJ2AccelerationModelActive )
        {
            acceleration = acceleration + astro::computeJ2Acceleration( gravitationalParameter,
                                                                        position,
                                                                        equatorialRadius,
                                                                        j2Coefficient );
        }

        // // Add solar radiation pressure acceleration if model is set to active.
        // if ( isRadiationPressureAccelerationModelActive )
        // {
        //     acceleration = sml::add( acceleration,
        //                              astro::computeJ2Acceleration( gravitationalParameter,
        //                                                            position,
        //                                                            equatorialRadius,
        //                                                            j2Coefficient ) );
        // }

        return State( { state[ 3 ],
                        state[ 4 ],
                        state[ 5 ],
                        acceleration[ 0 ],
                        acceleration[ 1 ],
                        acceleration[ 2 ] } );
    }

protected:

private:

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
};

} // namespace dustsim

#endif // DUSTSIM_DYNAMICAL_SYSTEM_HPP
