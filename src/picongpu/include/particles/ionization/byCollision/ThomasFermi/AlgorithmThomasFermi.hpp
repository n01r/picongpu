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
#include "particles/ionization/byCollision/ThomasFermi/TFFittingParameters.def"

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

            /* ionization condition */
            if (chargeState < protonNumber)
            {

                /* requires the "temperature" value in eV */
                float_X T_0 = temperature/math::pow(protonNumber,float_X(4./3.));

                float_X T_F = T_0 / (float_X(1.) + T_0);

                /* for all the fitting parameters @see TFFittingParameters.def */

                /** this is totally weird - I have to define temporary variables because
                 * otherwise the math::pow function won't recognize those at the
                 * exponent position */
                float_X TFA2_temp = TFA2;
                float_X TFA4_temp = TFA4;
                float_X TFBeta_temp = TFBeta;

                float_X   A = TFA1 * math::pow(T_0,TFA2_temp) + TFA3 * math::pow(T_0,TFA4_temp);

                float_X   B = -math::exp(TFB0 + TFB1*T_F + TFB2*math::pow(T_F,float_X(7.)));

                float_X   C = TFC1*T_F + TFC2;

                /* requires mass density in g/cm^3 */
                float_X   R = density/(protonNumber * A);

                float_X Q_1 = A * math::pow(R,B);

                float_X   Q = math::pow(math::pow(R,C) + math::pow(Q_1,C), float_X(1.)/C);

                float_X   x = TFAlpha * math::pow(Q,TFBeta_temp);

                /* Thomas-Fermi average ionization state */
                float_X ZStar = protonNumber * x / (float_X(1.) + x + math::sqrt(float_X(1.) + float_X(2.)*x));


                /* determine the new number of bound electrons from the TF ionization state */
                float_X newBoundElectrons = protonNumber - float_X(math::float2int_rn(ZStar));
                /* safety check to avoid double counting since recombination is not yet implemented */
                if (newBoundElectrons < parentIon[boundElectrons_])
                    /* update the particle attribute only if more free electrons are to be created */
                    parentIon[boundElectrons_] = newBoundElectrons;
            }

        }
    };

} // namespace ionization
} // namespace particles
} // namespace picongpu
