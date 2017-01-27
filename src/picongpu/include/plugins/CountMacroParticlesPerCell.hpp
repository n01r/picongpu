/**
 * Copyright 2017 Marco Garten
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

#include "mpi/MPIReduce.hpp"

namespace picongpu
{
using namespace PMacc;

template<class ParticlesType>
class CountMacroParticlesPerCell : public ILightweightPlugin
{

    ParticlesType *particles;
    std::ofstream outFile;
    /** @todo have only rank 0 create a file */
    bool writeToFile;
    std::string pluginName;
    std::string pluginPrefix;
    std::string filename;
    mpi::MPIReduce reduce;
    int mpiRank;

public:
    CountMacroParticlesPerCell() :
    pluginName("CountMacroParticles: count macro particles of a species in a cell"),
    pluginPrefix(ParticlesType::FrameType::getName() + std::string("_inCellCount")),
    filename(pluginPrefix + ".dat"),
    writeToFile(false)
    {
        /* register our plugin during creation */
        Environment<>::get().PluginConnector().registerPlugin(this);
    }

    std::string pluginGetName() const
    {
        return pluginName;
    }

    void notify(uint32_t currentStep)
    {
        countMacroParticles< CORE + BORDER >(currentStep);
    }

    void pluginRegisterHelp(po::options_description& desc)
    {
        /* register command line parameters for your plugin */
        desc.add_options()
            ((pluginPrefix+".period").c_str(), po::value<uint32_t > (&notifyPeriod)->default_value(0),
            "Enable CountMacroParticlesPerCell [for each n-th step]");
    }

    void setMappingDescription(MappingDesc *cellDescription)
    {
    }

private:
    uint32_t notifyPeriod;

    void pluginLoad()
    {
        if (notifyPeriod > 0)
        {
            /** determine rank 0 here to write output only once */
            writeToFile = reduce.hasResult(mpi::reduceMethods::Reduce());

            if (writeToFile)
            {
                outFile.open(filename.c_str(), std::ofstream::out | std::ostream::trunc);
                if (!outFile)
                {
                    std::cerr << "Can't open file [" << filename << "] for output, disable plugin output. " << std::endl;
                    writeToFile = false;
                }
                //create header of the file
                outFile << "#step count" << " \n";
            }
            /** Queue `notify()` of this plugin into the `notificationList` */
            Environment<>::get().PluginConnector().setNotificationPeriod(this, notifyPeriod);
        }
    }

    void pluginUnload()
    {
        /* called when plugin is unloaded, cleanup here */
        if (notifyPeriod > 0)
        {
            if (writeToFile)
            {
                outFile.flush();
                outFile << std::endl; //now all data are written to file
                if (outFile.fail())
                    std::cerr << "Error on flushing file [" << filename << "]. " << std::endl;
                outFile.close();
            }
        }
    }

    template< uint32_t AREA >
    void countMacroParticles(uint32_t currentStep)
    {
        /* calls a kernel that accesses the frames for the number of particles
         * in each cell, transfers that data to the host and writes it to a file
         */

        /** @todo call the kernel */
        /** @todo transfer data to host */
        /** @todo write data into file */
        if (writeToFile)
        {
            typedef std::numeric_limits< float_64 > dbl;

            /* dummy size */
            int rowSize = 8;

            outFile.precision(dbl::digits10);

            /* write data to file */
            outFile << currentStep << " "
                    << std::scientific; /*  for floating points, ignored for ints */

            for (int i = 0; i < rowSize; ++i)
            {
                outFile << std::scientific << double(1.0) << " ";
            }
            outFile << std::endl;
            /* endl: Flush any step to the file.
             * Thus, we will have data if the program should crash. */
        }
    }

};
} // namespace picongpu