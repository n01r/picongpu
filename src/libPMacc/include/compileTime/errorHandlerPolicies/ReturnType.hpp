/**
 * Copyright 2015-2017 Alexander Grund
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

#include "pmacc_types.hpp"

namespace PMacc
{
namespace errorHandlerPolicies
{

/** Returns the given type
 *  Binary meta function that takes any boost mpl sequence and a type
 */
template<typename T_ReturnType = bmpl::void_>
struct ReturnType
{
    template<typename T_MPLSeq, typename T_Value>
    struct apply
    {
        typedef T_ReturnType type;
    };
};

} // namespace errorHandlerPolicies
} // namespace PMacc
