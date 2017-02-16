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

namespace PMacc
{
struct KernelCountMacroParticlesPerCell
{
    template<class PBox, class T_DataBox, class Mapping>
    DINLINE void operator()(
        PBox pb,
        T_DataBox& inCellCounter,
        Mapping mapper
    ) const
    {

        typedef typename PBox::FramePtr FramePtr;
        const uint32_t Dim = Mapping::Dim;

        PMACC_SMEM( frame, FramePtr );
        PMACC_SMEM( counter, int );
        PMACC_SMEM( particlesInSuperCell, lcellId_t );


        typedef typename Mapping::SuperCellSize SuperCellSize;

        const DataSpace<Dim> threadIndex(threadIdx);
        const int linearThreadIdx = DataSpaceOperations<Dim>::template map<SuperCellSize > (threadIndex);
        const DataSpace<Dim> superCellIdx(mapper.getSuperCellIndex(DataSpace<Dim > (blockIdx)));

        /* multi-dim offset from the origin of the local domain on GPU
         * to the origin of the block of the in unit of cells
         */
        const DataSpace<Dim> blockCell = superCellIdx * SuperCellSize::toRT();

        /* subtract guarding cells to only have the simulation volume */
        const DataSpace<Dim> globalCellIndex = (superCellIdx * SuperCellSize::toRT() + threadIndex) - mapper.getGuardingSuperCells() * SuperCellSize::toRT();

        inCellCounter(globalCellIndex) = pb.getCellCount(superCellIdx)[linearThreadIdx];

//        if (linearThreadIdx == 0)
//        {
//            frame = pb.getLastFrame(superCellIdx);
//            particlesInSuperCell = pb.getSuperCell(superCellIdx).getSizeLastFrame();
//            counter = 0;
//        }
//        __syncthreads();
//        if (!frame.isValid())
//            return; //end kernel if we have no frames
//        filter.setSuperCellPosition((superCellIdx - mapper.getGuardingSuperCells()) * mapper.getSuperCellSize());
//        while (frame.isValid())
//        {
//            if (linearThreadIdx < particlesInSuperCell)
//            {
//                if (filter(*frame, linearThreadIdx))
//                    nvidia::atomicAllInc(&counter);
//            }
//            __syncthreads();
//            if (linearThreadIdx == 0)
//            {
//                frame = pb.getPreviousFrame(frame);
//                particlesInSuperCell = math::CT::volume<SuperCellSize>::type::value;
//            }
//            __syncthreads();
//        }

        __syncthreads();
//        if (linearThreadIdx == 0)
//        {
//            atomicAdd(gCounter, (uint64_cu) counter);
//        }
    }
};

struct NumberMacroParticlesPerCell
{

    /** Get particle count
     *
     * @tparam AREA area were particles are counted (CORE, BORDER, GUARD)
     *
     * @param buffer source particle buffer
     * @param cellDescription instance of MappingDesction
     * @param filter filter instance which must inherit from PositionFilter
     * @return number of particles in defined area
     */
    template<uint32_t AREA, class PBuffer, class CellDesc, class T_inCellCounter>
    static void countOnDevice(PBuffer& buffer, CellDesc cellDescription, T_inCellCounter& inCellCounter)
    {

        auto block = CellDesc::SuperCellSize::toRT();

        AreaMapping<AREA, CellDesc> mapper(cellDescription);

        printf("Pointer address of inCellCounter #3 %p : \n", &inCellCounter);

        PMACC_KERNEL(KernelCountMacroParticlesPerCell{})
            (mapper.getGridDim(), block)
            (buffer.getDeviceParticlesBox(),
             inCellCounter.getDeviceBuffer().getDataBox(),
             mapper);

        inCellCounter.deviceToHost();
    }
};

} // namespace PMacc

namespace picongpu
{
using namespace PMacc;

template<class ParticlesType>
class CountMacroParticlesPerCell : public ILightweightPlugin
{

private:
    uint32_t notifyPeriod;
    GridBuffer<uint32_t, simDim>* inCellCounter;
    ParticlesType *particles;
    MappingDesc *cellDescription;
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
    writeToFile(false),
    particles(nullptr),
    cellDescription(nullptr),
    inCellCounter(nullptr)
    {
        /* register our plugin during creation */
        Environment<>::get().PluginConnector().registerPlugin(this);
    }

    virtual ~CountMacroParticlesPerCell()
    {

    }

    std::string pluginGetName() const
    {
        return pluginName;
    }

    void notify(uint32_t currentStep)
    {
        DataConnector &dc = Environment<>::get().DataConnector();

        particles = &(dc.getData< ParticlesType > (ParticlesType::FrameType::getName(), true));

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
        this->cellDescription = cellDescription;
    }

    void pluginLoad()
    {
        if (notifyPeriod > 0)
        {
            /** determine rank 0 here to write output only once */
            writeToFile = reduce.hasResult(mpi::reduceMethods::Reduce());

            inCellCounter = new GridBuffer<uint32_t, simDim>(cellDescription->getGridLayout());
            printf("Pointer address of inCellCounter #1 %p : \n", inCellCounter);

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
            __delete(inCellCounter);
        }
    }

    template< uint32_t AREA >
    void countMacroParticles(uint32_t currentStep)
    {
        printf("Pointer address of inCellCounter #2 %p : \n", inCellCounter);
        inCellCounter->getHostBuffer().setValue(0);
        inCellCounter->getDeviceBuffer().setValue(0);

        /*count local particles*/
        PMacc::NumberMacroParticlesPerCell::countOnDevice<AREA>(
            *particles,
            *cellDescription,
            *inCellCounter
        );

        unsigned int DataSpace<SIMDIM> myIndex;
        myIndex[0] = 1;
        myIndex[1] = 1;
        myIndex[2] = 1;

        printf("inCellCounter Value #4 %u : \n", inCellCounter(myIndex));

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
