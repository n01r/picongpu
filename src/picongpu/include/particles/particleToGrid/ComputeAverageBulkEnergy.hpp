/**
 * Copyright 2013-2016 Axel Huebl, Heiko Burau, Rene Widera
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

#include "simulation_defines.hpp"
#include "pmacc_types.hpp"

#include "math/Vector.hpp"
#include "particles/particleToGrid/ComputeGridValuePerFrame.def"
#include "particles/particleToGrid/derivedAttributes/DerivedAttributes.hpp"

#include "algorithms/Gamma.hpp"

#include <vector>
#include "cuSTL/container/DeviceBuffer.hpp"

namespace picongpu
{
namespace particleToGrid
{

    template<typename T_Species>
    struct ComputeAverageBulkEnergy
    {


        typedef typename CreateEnergyOperation<T_Species>::type::Solver EnergySolver;
        typedef typename CreateCounterOperation<T_Species>::type::Solver CounterSolver;

        container::DeviceBuffer<float_X,SIMDIM>* counterBuffer;

        const SubGrid<simDim>& subGrid = Environment<simDim>::get().SubGrid();
        DataSpace<simDim> localDomain = subGrid.getLocalDomain( ).size;

        HDINLINE ComputeAverageBulkEnergy()
        {
            counterBuffer = new container::DeviceBuffer<float_X,SIMDIM>(localDomain);
        }



        DataConnector &dc = Environment<>::get().DataConnector();
        /*load FieldTmp without copy data to host*/
        FieldTmp* fieldTmp = &(dc.getData<FieldTmp > (FieldTmp::getName(), true));

        /* load species without copying the particle data to the host */
        T_Species* speciesTmp = &(dc.getData<T_Species >(T_Species::FrameType::getName(), true));

        ~ComputeAverageBulkEnergy()
        {
            dc.releaseData(T_Species::FrameType::getName());
        }

        template<typename T_Species, typename FieldTmp>
        DINLINE void operator()(T_Species* speciesTmp, FieldTmp* fieldTmp, counterBuffer, uint32_t currentStep)
        {

            const uint32_t area = bmpl::int_<CORE+BORDER>;

            fieldTmp->getGridBuffer().getDeviceBuffer().setValue(FieldTmp::ValueType(0.0));
            /* count the contributing particles in each cell */
            fieldTmp->computeValue < area, CounterSolver > (*speciesTmp, currentStep);

            /* cast libPMacc Buffer to cuSTL Buffer and subtract guarding cells */
            BOOST_AUTO(fieldTmp_coreBorder,
                         fieldTmp->getGridBuffer().
                         getDeviceBuffer().cartBuffer().
                         view(this->cellDescription->getGuardingSuperCells()*BlockDim::toRT(),
                              this->cellDescription->getGuardingSuperCells()*-BlockDim::toRT()));

            /* temporarily save the number of the macroparticles for later averaging */
            *counterBuffer = fieldTmp_coreBorder;

            fieldTmp->getGridBuffer().getDeviceBuffer().setValue(FieldTmp::ValueType(0.0));
            /* sum up the weighted particle energies in each cell */
            fieldTmp->computeValue < area, EnergySolver > (*speciesTmp, currentStep);

            /* cast libPMacc Buffer to cuSTL Buffer and subtract guarding cells */
            BOOST_AUTO(fieldTmp_coreBorder,
                         fieldTmp->getGridBuffer().
                         getDeviceBuffer().cartBuffer().
                         view(this->cellDescription->getGuardingSuperCells()*BlockDim::toRT(),
                              this->cellDescription->getGuardingSuperCells()*-BlockDim::toRT()));



            using namespace lambda;
            using namespace PMacc::math::math_functor;

            algorithm::kernel::Foreach<BlockDim>()(fieldTmp_coreBorder.zone(), fieldTmp_coreBorder.origin(),*counterBuffer.origin(),__1 = __1 / __2)

        }
    };

} // namespace particleToGrid
} // namespace picongpu

