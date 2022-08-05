/*
 * Copyright (c) 2009-2022 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#pragma once

#include <iomanip>
#include <iostream>
#include <string>

#include "dustsim/typedefs.hpp"

namespace dustsim
{

//! Print value to stream.
/*!
 * Prints a specified value to stream provided, given a specified width and a filler character.
 *
 * @tparam      DataType  Type for specified value
 * @param[out]  stream    Output stream
 * @param[in]   value     Specified value to print
 * @param[in]   width     Width of value printed to stream, in terms of number of characters
 *                       (default = 25)
 * @param[in]   filler    Character used to fill fixed-width (default = ' ')
 */
template <typename DataType>
inline void print(std::ostream& stream,
                  const DataType value,
                  const int width = 25,
                  const char filler = ' ')
{
    stream << std::left << std::setw(width) << std::setfill(filler) << value;
}

//! Print metadata parameter to stream.
/*!
 * Prints metadata parameter to stream provided, given a specified name, value, units, and
 * delimiter.
 *
 * @tparam      DataType       Type for specified value
 * @param[out]  stream         Output stream
 * @param[in]   parameterName  Name for metadata parameter
 * @param[in]   value          Specified value to print
 * @param[in]   units          Units for value
 * @param[in]   delimiter      Character used to delimiter entries in stream (default = ' ')
 * @param[in]   width          Width of value printed to stream, in terms of number of characters
 *                             (default = 25)
 * @param[in]  filler          Character used to fill fixed-width (default = ' ')
 */
template <typename DataType>
inline void print(std::ostream& stream,
                  const std::string& parameterName,
                  const DataType value,
                  const std::string& units,
                  const char delimiter = ',',
                  const int width = 25,
                  const char filler = ' ')
{
    print(stream, parameterName, width, filler);
    stream << delimiter;
    print(stream, value, width, filler);
    stream << delimiter;
    print(stream, units, width, filler);
}

} // namespace dustsim
