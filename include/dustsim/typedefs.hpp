/*
 * Copyright (c) 2009-2025 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#pragma once

#include <vector>

namespace dustsim
{

//! Set type for integer numbers.
typedef int Int;

//! Set type for floating-point real numbers.
typedef double Real;

//! Set type for n-dimensional vector.
typedef std::vector<Real> Vector;

//! Define numerical integrators.
enum Integrator
{
  rk4,
  rkf45,
  rkf78
};

} // namespace dustsim
