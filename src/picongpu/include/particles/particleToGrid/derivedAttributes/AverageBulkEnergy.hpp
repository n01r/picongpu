/**
 * Copyright 2016 Marco Garten, Axel Huebl, Heiko Burau, Rene Widera,
 *                Felix Schmitt
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

#include "particles/particleToGrid/derivedAttributes/AverageBulkEnergy.def"

#include "simulation_defines.hpp"



namespace picongpu
{
namespace particles
{
namespace ionization
{


/* functor to calculate average bulk energy */
template<typename T_SpeciesName, typename T_Area>
struct ComputeAverageBulkEnergy
{
    typedef typename T_SpeciesName::type SpeciesName;
    static const uint32_t area = T_Area::value;

    HINLINE void operator()( FieldTmp* fieldTmp,
                             FieldTmp* counterBuffer,
                             const uint32_t currentStep) const
    {
        DataConnector &dc = Environment<>::get().DataConnector();

        /* load species without copying the particle data to the host */
        SpeciesName* speciesTmp = &(dc.getData<SpeciesName >(SpeciesName::FrameType::getName(), true));

        /* run algorithm */
        typedef typename particleToGrid::CreateAveragedEnergyOperation<SpeciesName>::type::Solver AverageBulkEnergySolver;
        fieldTmp->computeValue < area, AverageBulkEnergySolver > (*speciesTmp, currentStep);




        Hier weiter!
        jehfiluheuöahguöahrgöharghaörghaörggABlärhx!



        

        dc.releaseData(SpeciesName::FrameType::getName());
    }
};

typedef SuperCellSize BlockDim;

DataConnector &dc = Environment<>::get().DataConnector();

/* load FieldTmp without copying data to host */
FieldTmp* fieldTmp = &(dc.getData<FieldTmp > (FieldTmp::getName(), true));
/* reset average bulk energy values to zero */
fieldTmp->getGridBuffer().getDeviceBuffer().setValue(FieldTmp::ValueType(0.0));

/* calculate and add the average bulk energy from the electron species in FieldTmp */
/* @todo filter for species and make the filter parameters user-customizable */
ForEach<Species1, picongpu::particles::ionization::ComputeAverageBulkEnergy<bmpl::_1,bmpl::int_<CORE + BORDER> >, MakeIdentifier<bmpl::_1> > computeAverageBulkEnergy;
computeAverageBulkEnergy(forward(fieldTmp), currentStep);

/* add results of all species that are still in GUARD to next GPUs BORDER */
EventTask fieldTmpEvent = fieldTmp->asyncCommunication(__getTransactionEvent());
__setTransactionEvent(fieldTmpEvent);

} //namespace ionization
} //namespace particles
} // namespace picongpu
