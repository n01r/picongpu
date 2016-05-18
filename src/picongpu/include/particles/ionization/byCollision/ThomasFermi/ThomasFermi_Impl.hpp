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
#include "math/vector/Size_t.hpp"
#include "simulation_defines.hpp"
#include "traits/Resolve.hpp"
#include "mappings/kernel/AreaMapping.hpp"

#include "fields/FieldTmp.hpp"

#include "particles/ionization/byCollision/ThomasFermi/ThomasFermi.def"
#include "particles/ionization/byCollision/ThomasFermi/AlgorithmThomasFermi.hpp"
#include "particles/ionization/ionization.hpp"

#include "compileTime/conversion/TypeToPointerPair.hpp"
#include "memory/boxes/DataBox.hpp"

#include "particles/ParticlesFunctors.hpp"

#include "particles/particleToGrid/ComputeGridValuePerFrame.def"
#include "particles/particleToGrid/ComputeGridValuePerFrame.def"

namespace picongpu
{
namespace particles
{
namespace ionization
{

    /** \struct ThomasFermi_Impl
     *
     * \brief Thomas-Fermi pressure ionization - Implementation
     *
     * \tparam T_DestSpecies electron species to be created
     * \tparam T_SrcSpecies particle species that is ionized
     */
    template<typename T_DestSpecies, typename T_SrcSpecies>
    struct ThomasFermi_Impl
    {

        typedef T_DestSpecies DestSpecies;
        typedef T_SrcSpecies  SrcSpecies;

        typedef typename SrcSpecies::FrameType FrameType;

        /** Solvers for electron macroparticle energy and electron macroparticle weighting
         *
         * Those are used to calculate the average bulk electron energy
         * with respect to the cell the ion is located in
         *
         * \todo Get the bulk energy of all electrons in the simulation,
         *       not only the of the destination species
         */
        //typedef typename particleToGrid::CreateEnergyOperation<T_DestSpecies>::type::Solver EnergySolver;
        //typedef typename particleToGrid::CreateWeightingOperation<T_DestSpecies>::type::Solver WeightingSolver;

        /** Solver for density of the ion species
         *
         *  \todo Include all ion species because the model requires the
         *        density of ionic potential wells
         *  \todo Implement the real density and rename density to charge density
         *        because this is still what is implemented now
         */
        //typedef typename particleToGrid::CreateDensityOperation<T_SrcSpecies>::type::Solver DensitySolver;

        private:

            /* define ionization ALGORITHM (calculation) for ionization MODEL */
            typedef particles::ionization::AlgorithmThomasFermi IonizationAlgorithm;

            typedef MappingDesc::SuperCellSize TVec;
//
//            /* "temperature" value type */
//            typedef FieldTmp::ValueType ValueType_T;
//            /* global memory "temperature"-field device databox */
//            FieldTmp::DataBoxType temperatureBox;
//
//            /* shared memory "temperature"-field device databox */
//            PMACC_ALIGN(cachedT, DataBox<SharedBox<ValueType_T, typename BlockArea::FullSuperCellSize,0> >);

        public:
            /* host constructor */
            ThomasFermi_Impl(const uint32_t currentStep)
            {
//                DataConnector &dc = Environment<>::get().DataConnector();
//                /* initialize pointers on host-side "temperature"-field databox */
//                FieldTmp* fieldTmp = &(dc.getData<FieldTmp > (FieldTmp::getName(), true));
//                /* initialize device-side "temperature"-field databox */
//                temperatureBox = fieldTmp->getDeviceDataBox();
//
//                AbbbEhoiwehfoiwhefilhnalihförngösdfngsöldfnbsöjdfbnadkfnasködfbnskdfbnskdf.bnsd.kbnsdk.gbnsk.dgbnskd.gbn
//                jdsfnvjsdhnfuiovgahruglhfuivhalfiuvnaifnviafnviefhnvaihvnafhnvafnva
//                Hier weitermachen!

            }

            /** Initialization function on device
             *
             * \brief Cache EM-fields on device
             *         and initialize possible prerequisites for ionization, like e.g. random number generator.
             *
             * This function will be called inline on the device which must happen BEFORE threads diverge
             * during loop execution. The reason for this is the `__syncthreads()` call which is necessary after
             * initializing the E-/B-field shared boxes in shared memory.
             *
             * @param blockCell Offset of the cell from the origin of the local domain
             *                  <b>including guarding supercells</b> in units of cells
             * @param linearThreadIdx Linearized thread ID inside the block
             * @param localCellOffset Offset of the cell from the origin of the local
             *                        domain, i.e. from the @see BORDER
             *                        <b>without guarding supercells</b>
             */
            DINLINE void init(const DataSpace<simDim>& blockCell, const int& linearThreadIdx, const DataSpace<simDim>& localCellOffset)
            {


            }

            /** Functor implementation
             *
             * \param ionFrame reference to frame of the to-be-ionized particles
             * \param localIdx local (linear) index in super cell / frame
             * \param newMacroElectrons reference to variable for each thread that stores the number
             *        of macro electrons to be created during the current time step
             */
            DINLINE void operator()(FrameType& ionFrame, int localIdx, unsigned int& newMacroElectrons)
            {
                /* alias for the single macro-particle */
                PMACC_AUTO(particle,ionFrame[localIdx]);
                /* particle position, used for field-to-particle interpolation */
                floatD_X pos = particle[position_];
                const int particleCellIdx = particle[localCellIdx_];
                /* multi-dim coordinate of the local cell inside the super cell */
                DataSpace<TVec::dim> localCell(DataSpaceOperations<TVec::dim>::template map<TVec > (particleCellIdx));

                /* define number of bound macro electrons before ionization */
                float_X prevBoundElectrons = particle[boundElectrons_];

                /* density dummy for testing / g/cm^3*/
                float_X density = 10.0;
                /* "temperature" dummy for testing / eV */
                float_X temperature = 10.0;

                /* this is the point where actual ionization takes place */
                IonizationAlgorithm ionizeAlgo;
                ionizeAlgo(
                    temperature, density,
                    particle
                    );

                /* determine number of new macro electrons to be created */
                newMacroElectrons = prevBoundElectrons - particle[boundElectrons_];

            }

    };

} // namespace ionization
} // namespace particles
} // namespace picongpu
