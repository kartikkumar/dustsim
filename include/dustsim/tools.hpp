/*
 * Copyright (c) 2017, K. Kumar, Delft University of Technology (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef DUSTSIM_TOOLS_HPP
#define DUSTSIM_TOOLS_HPP

#include <iostream>
#include <string>

#include <rapidjson/document.h>

#include <astro/astro.hpp>

#include "dustsim/typedefs.hpp"

namespace dustsim
{

//! Find parameter.
/*!
 * Finds parameter in config stored in JSON document. An error is thrown if the parameter cannot
 * be found. If the parameter is found, an iterator to the member in the JSON document is returned.
 *
 * @param[in] config        JSON document containing config parameters
 * @param[in] parameterName Name of parameter to find
 * @return                  Iterator to parameter retreived from JSON document
 */
ConfigIterator find( const rapidjson::Document& config, const std::string& parameterName );

//! Class to write state history to an output stream.
/*!
 * This class contains parameters and functions to write the state history to a user-provided
 * output stream, e.g., console, file.
 */
class StateHistoryWriter
{
public:

    //! Construct state history writer.
    /*!
     * Constructure state history writer to write state history to a user-defined output stream.
     *
     * @param[in] aStateHistoryStream   Output stream
     */

    StateHistoryWriter( std::ostream& aStateHistoryStream,
                        const Real aGravitationalParameter )
        : stateHistoryStream( aStateHistoryStream ),
          gravitationalParameter( aGravitationalParameter )
    { }

    //! Overload ()-operator to write state to output stream.
    /*!
     * Overloads ()-operator to write current state state and epoch to a user-defined output stream.
     * The output stream is defined through the class constructor.
     *
     * This function fulfills the prototype required to work with the Boost Odeint library.
     *
     * @sa boost::odeint:integrator
     * @param[in] state Current state
     * @param[in] time  Current epoch
     */
    void operator( )( const State& state , const double time )
    {
        const State stateInKeplerElements
            = astro::convertCartesianToKeplerianElements( state, gravitationalParameter );

        stateHistoryStream << time << ','
                           << state[ 0 ] << ','
                           << state[ 1 ] << ','
                           << state[ 2 ] << ','
                           << state[ 3 ] << ','
                           << state[ 4 ] << ','
                           << state[ 5 ] << ','
                           << stateInKeplerElements[ 0 ] << ','
                           << stateInKeplerElements[ 1 ] << ','
                           << stateInKeplerElements[ 2 ] << ','
                           << stateInKeplerElements[ 3 ] << ','
                           << stateInKeplerElements[ 4 ] << ','
                           << stateInKeplerElements[ 5 ]
                           << std::endl;
    }

protected:

private:

    //! Output stream to write state history to.
    std::ostream& stateHistoryStream;

    //! Gravitation parameter of central body [km^3 s^-2].
    const Real gravitationalParameter;
};

} // namespace dustsim

#endif // DUSTSIM_TOOLS_HPP
