/*
 * Copyright (c) 2009-2018, K. Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef DUSTSIM_STATE_HPP
#define DUSTSIM_STATE_HPP

#include <iostream>
#include <map>

#include "dustsim/typedefs.hpp"

namespace dustsim
{

class State
{
public:

    State( const Vector& aVector )
        : vector( aVector )
    { }

    const int size( ) const { return vector.size( ); }
    const Real operator[ ] ( const int i ) const { return vector[ i ]; }
    Real& operator[ ] ( const int i ) { return vector[ i ]; }

    State& operator=( const State& rightHandSide );
    State& operator+=( const State& rightHandSide );

    bool operator==( const State& rightHandSide ) const;
    bool operator<( const State& rightHandSide ) const;
    bool operator<( const Real rightHandSide ) const;
    bool operator>( const State& rightHandSide ) const;
    bool operator>( const Real rightHandSide ) const;

    friend State operator+( const State& leftHandSide, const State& rightHandSide );
    friend State operator*( const Real multiplier, const State& state );
    friend State operator*( const State& state, const Real multiplier );

    friend std::ostream& operator<<( std::ostream& stream, const State& state );

protected:
private:

    Vector vector;
};

typedef std::map< Real, State > StateHistory;

} // namespace dustsim

#endif // DUSTSIM_STATE_HPP
