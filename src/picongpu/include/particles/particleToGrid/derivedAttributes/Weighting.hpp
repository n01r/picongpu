/**
 * Copyright 2013-2016 Axel Huebl, Rene Widera, Marco Garten
 *
 * This file is part of PIConGPU.
 *
 * PIConGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PIConGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PIConGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "particles/particleToGrid/derivedAttributes/ChargeDensity.def"

#include "simulation_defines.hpp"

namespace picongpu
{
namespace particleToGrid
{
namespace derivedAttributes
{

    HDINLINE float1_64
    Weighting::getUnit() const
    {
        return 0.0;
    }

    template< class T_Particle >
    DINLINE float_X
    Weighting::operator()( T_Particle& particle ) const
    {
        /* read existing attributes */
        const float_X weighting = particle[weighting_];

        /* calculate new attribute */
        const float_X particleWeighting = weighting;

        /** return attribute particleWeighting
         *
         * Even though there already is an attribute "weighting" it is required
         * to have a "weighted weighting", i.e. weighting*AssignmentFunction to
         * account for the macroparticle shape. \see ComputeGridValuePerFrame
         * This is needed to calculate an
         * average bulk energy ("temperature") for the Thomas-Fermi ionization
         * model.
         */
        return particleWeighting;
    }
} /* namespace derivedAttributes */
} /* namespace particleToGrid */
} /* namespace picongpu */
