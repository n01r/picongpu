/**
 * Copyright 2016 Marco Garten
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PIConGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "pmacc_types.hpp"
#include "simulation_defines.hpp"
#include "particles/traits/GetAtomicNumbers.hpp"
#include "traits/attribute/GetChargeState.hpp"
#include "algorithms/math/floatMath/floatingPoint.tpp"

/** \file AlgorithmThomasFermi.hpp
 *
 * IONIZATION ALGORITHM for the Thomas-Fermi model
 *
 * - implements the calculation of average ionization degree and changes charge states
 *   by decreasing the number of bound electrons
 * - is called with the IONIZATION MODEL, specifically by setting the flag in @see speciesDefinition.param
 */

namespace picongpu
{
namespace particles
{
namespace ionization
{

    /** \struct AlgorithmThomasFermi
     *
     * \brief calculation for the Thomas-Fermi pressure ionization model
     *
     * This model uses local density and "temperature" values as input
     * parameters. To be able to speak of a "temperature" an equilibrium state
     * would be required. Typical high power laser-plasma interaction is highly
     * non-equilibrated, though. The name "temperature" is kept to illustrate
     * the origination from the Thomas-Fermi model. It is nevertheless
     * more accurate to think of it as an averaged kinetic energy.
     */
    struct AlgorithmThomasFermi
    {

        /** Functor implementation
         *
         * \tparam DensityType type of mass density
         * \tparam TemperatureType type of temperature
         * \tparam ParticleType type of particle to be ionized
         *
         * \param density mass density value at t=0
         * \param temperature temperature value at t=0
         * \param parentIon particle instance to be ionized with position at t=0 and momentum at t=-1/2
         */
        template<typename DensityType, typename TemperatureType, typename ParticleType >
        HDINLINE void
        operator()( const TemperatureType temperature, const DensityType density, ParticleType& parentIon )
        {

            const float_X protonNumber = GetAtomicNumbers<ParticleType>::type::numberOfProtons;
            float_X chargeState = attribute::getChargeState(parentIon);

            uint32_t cs = math::float2int_rd(chargeState);

            /* ionization condition */
            if (math::abs(density) / ATOMIC_UNIT_EFIELD >= AU::IONIZATION_EFIELD_HYDROGEN[cs] && chargeState < protonNumber)
            {
                /* set new particle charge state */
                parentIon[boundElectrons_] -= float_X(1.0);
            }

        }
    };

} // namespace ionization
} // namespace particles
} // namespace picongpu
