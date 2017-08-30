/*
 * Copyright (c) 2017, K. Kumar, Delft University of Technology (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef DUSTSIM_TYPEDEFS_HPP
#define DUSTSIM_TYPEDEFS_HPP

#include <boost/array.hpp>

#include <rapidjson/document.h>

namespace dustsim
{

//! Set type for floating-point real numbers.
typedef double Real;

//! Set type for 3-dimensional vector.
typedef boost::array< Real, 3 > Vector3;

//! Set type for 6-dimensional vector.
typedef boost::array< Real, 6 > Vector6;

//! Set type for position vector to 3D vector.
typedef Vector3 Position;

//! Set type for velocity vector to 3D vector.
typedef Vector3 Velocity;

//! Set type for state vector to 6D vector.
typedef Vector6 State;

//! JSON config iterator.
typedef rapidjson::Value::ConstMemberIterator ConfigIterator;

} // namespace dustsim

#endif // DUSTSIM_TYPEDEFS_HPP