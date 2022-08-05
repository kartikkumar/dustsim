/*
 * Copyright (c) 2009-2022 Kartik Kumar (me@kartikkumar.com)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <iostream>

#include "dustsim/state.hpp"

namespace dustsim
{

State& State::operator=(const State& rightHandSide)
{
    // Check for self-assignment.
    if (this != &rightHandSide)
    {
        vector = rightHandSide.vector;
    }

    return *this;
}

State& State::operator+=(const State& rightHandSide)
{
    *this = *this + rightHandSide;
    return *this;
}

bool State::operator==(const State& rightHandSide) const
{
    return ((*this).vector == rightHandSide.vector);
}

bool State::operator<(const State& rightHandSide) const
{
    return ((*this).vector < rightHandSide.vector);
}

bool State::operator<(const Real rightHandSide) const
{
    bool isLess = true;
    for (int i = 0; i < (*this).size(); ++i)
    {
        if ((*this).vector[i] > rightHandSide) {isLess = false;};
        break;
    }
    return isLess;
}

bool State::operator>(const State& rightHandSide) const
{
    return ((*this).vector > rightHandSide.vector);
}

bool State::operator>(const Real rightHandSide) const
{
    bool isMore = true;
    for (int i = 0; i < (*this).size(); ++i)
    {
        if ((*this).vector[i] < rightHandSide) {isMore = false;};
        break;
    }
    return isMore;
}

State operator+(const State& leftHandSide, const State& rightHandSide)
{
    Vector vector(leftHandSide.size());
    for (unsigned int i = 0; i < vector.size(); i++)
    {
        vector[i] = leftHandSide[i] + rightHandSide[i];
    }

    return State(vector);
}

State operator*(const Real multiplier, const State& state)
{
    Vector vector(state.size());
    for (unsigned int i = 0; i < vector.size(); i++)
    {
        vector[i] = multiplier * state[i];
    }
    return State(vector);
}

State operator*(const State& state, const Real multiplier)
{
    return multiplier * state;
}

std::ostream& operator<<(std::ostream& stream, const State& state)
{
    for (int i = 0; i < state.size() - 1; i++) {stream << state[i] << ",";}
    stream << state[state.size() - 1];
    return stream;
}

} // namespace dustsim
