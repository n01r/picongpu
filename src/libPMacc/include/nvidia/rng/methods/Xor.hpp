/**
 * Copyright 2013-2017 Heiko Burau, Rene Widera
 *
 * This file is part of libPMacc.
 *
 * libPMacc is free software: you can redistribute it and/or modify
 * it under the terms of either the GNU General Public License or
 * the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libPMacc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License and the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * and the GNU Lesser General Public License along with libPMacc.
 * If not, see <http://www.gnu.org/licenses/>.
 */


#pragma once

#include <curand_kernel.h>
#include "pmacc_types.hpp"

namespace PMacc
{
namespace nvidia
{
namespace rng
{
namespace methods
{

class Xor
{
public:
    typedef curandStateXORWOW_t StateType;
    typedef StateType* StatePtr;

    HDINLINE Xor()
    {
    }

    DINLINE Xor(uint32_t seed, uint32_t subsequence = 0, uint32_t offset = 0)
    {
        curand_init(seed, subsequence, offset, &state);
    }

    HDINLINE Xor(const Xor& other): state(other.state)
    {

    }

protected:

    DINLINE curandStateXORWOW_t* getStatePtr()
    {
        return &state;
    }

    DINLINE curandStateXORWOW_t& getState()
    {
        return state;
    }

private:
    PMACC_ALIGN(state, StateType);
};
}
}
}
}
