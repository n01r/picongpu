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

namespace picongpu
{
using namespace PMacc;

template<class T_SpeciesType>
class AverageBulkEnergy : public ISimulationPlugin
{
    private:
        typedef T_SpeciesType SpeciesType;

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

        void notify(uint32_t currentStep)
        {
            /* notification callback for simulation step currentStep
             * called every notifyPeriod steps */
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
