/**
 * Copyright 2013-2016 Marco Garten, Axel Huebl, Felix Schmitt, Heiko Burau,
 *                     Rene Widera, Richard Pausch, Benjamin Worpitz
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

#include "plugins/ILightweightPlugin.hpp"
#include "dataManagement/DataConnector.hpp"
#include "fields/FieldTmp.hpp"

namespace picongpu
{
using namespace PMacc;

namespace detail
{
    // functor for all species to calculate density
    template<typename T_SpeciesName, typename T_Area>
    struct ComputeAverageBulkEnergy
    {
        typedef typename T_SpeciesName::type SpeciesName;
        static const uint32_t area = T_Area::value;

        HINLINE void operator()( FieldTmp* fieldTmp,
                                 const uint32_t currentStep) const
        {
            DataConnector &dc = Environment<>::get().DataConnector();

            /* load species without copying the particle data to the host */
            SpeciesName* speciesTmp = &(dc.getData<SpeciesName >(SpeciesName::FrameType::getName(), true));

            /* run algorithm */
            typedef typename CreateEnergyOperation<SpeciesName>::type::Solver EnergySolver;
            fieldTmp->computeValue < area, EnergySolver > (*speciesTmp, currentStep);
            dc.releaseData(SpeciesName::FrameType::getName());
        }
    };
} // namespace detail

template<class T_SpeciesType>
class AverageBulkEnergy : public ISimulationPlugin
{
    private:
        typedef T_SpeciesType SpeciesType;

        /* only rank 0 creates a file */
        bool writeToFile;

    public:
        AverageBulkEnergy()
        {
            /* register our plugin during creation */
            Environment<>::get().PluginConnector().registerPlugin(this);
        }

        std::string pluginGetName() const
        {
            return "AverageBulkEnergy";
        }

        /* notification callback for simulation step currentStep
         * called every notifyPeriod steps */
        void notify(uint32_t currentStep)
        {
            typedef SuperCellSize BlockDim;

            DataConnector &dc = Environment<>::get().DataConnector();

            /* load FieldTmp without copy data to host */
            FieldTmp* fieldTmp = &(dc.getData<FieldTmp > (FieldTmp::getName(), true));
            /* reset density values to zero */
            fieldTmp->getGridBuffer().getDeviceBuffer().setValue(FieldTmp::ValueType(0.0));

            /* calculate the particle energies for all particles */
            ForEach<VectorAllSpecies, picongpu::detail::ComputeAverageBulkEnergy<bmpl::_1,bmpl::int_<CORE + BORDER> >, MakeIdentifier<bmpl::_1> > computeAverageBulkEnergy;
            computeAverageBulkEnergy(forward(fieldTmp), currentStep);

            /* add results of all species that are still in GUARD to next GPUs BORDER */
            EventTask fieldTmpEvent = fieldTmp->asyncCommunication(__getTransactionEvent());
            __setTransactionEvent(fieldTmpEvent);
            std::cout << "Hello World!" << std::endl;

        }

        void pluginRegisterHelp(po::options_description& desc)
        {
            /* register command line parameters for your plugin */
            desc.add_options()
              ("avgBulkEnergy.period", po::value<uint32_t > (&notifyPeriod)->default_value(0),
               "Enable AverageBulkEnergy [for each n-th step]");
        }

        void setMappingDescription(MappingDesc *cellDescription)
        {
        }

        void restart(uint32_t restartStep, const std::string restartDirectory)
        {
            /* restart from a checkpoint here
             * will be called only once per simulation and before notify() */
        }

        void checkpoint(uint32_t currentStep, const std::string restartDirectory)
        {
        /* create a persistent checkpoint here
         * will be called before notify() if both will be called for the same timestep */
        }

    private:
        uint32_t notifyPeriod;

        void pluginLoad()
        {
            /* called when plugin is loaded, command line flags are available here
             * set notification period for our plugin at the PluginConnector */
            Environment<>::get().PluginConnector().setNotificationPeriod(this, notifyPeriod);
        }

        void pluginUnload()
        {
            /* called when plugin is unloaded, cleanup here */
        }
    };
} // namespace picongpu
