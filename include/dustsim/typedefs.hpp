/*
 * Copyright (c) 2009-2018, K. Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#ifndef DUSTSIM_TYPEDEFS_HPP
#define DUSTSIM_TYPEDEFS_HPP

#include <vector>

#include <rapidjson/document.h>

namespace dustsim
{

//! Set type for integer numbers.
typedef int Int;

//! Set type for floating-point real numbers.
typedef double Real;

//! Set type for n-dimensional vector.
typedef std::vector< Real > Vector;

//! JSON config iterator.
typedef rapidjson::Value::ConstMemberIterator ConfigIterator;

//! Define numerical integrators.
enum Integrator
{
  rk4,
  rkf45,
  rkf78
};

} // namespace dustsim

#endif // DUSTSIM_TYPEDEFS_HPP
