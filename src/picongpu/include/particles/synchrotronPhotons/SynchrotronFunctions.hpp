/**
 * Copyright 2015-2017 Heiko Burau
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

#include "simulation_defines.hpp"

#include "cuSTL/container/HostBuffer.hpp"
#include "cuSTL/cursor/Cursor.hpp"
#include "cuSTL/cursor/navigator/PlusNavigator.hpp"
#include "cuSTL/cursor/tools/LinearInterp.hpp"
#include "cuSTL/cursor/BufferCursor.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/math/tr1.hpp> /* cyl_bessel_k */

namespace picongpu
{
namespace particles
{
namespace synchrotronPhotons
{

namespace detail
{

/** Map `x` to the internal lookup table and return the result of the
 * first or the second synchrotron function for `x`.
 */
struct MapToLookupTable
{
    typedef typename ::PMacc::result_of::Functor<
        ::PMacc::cursor::tools::LinearInterp<float_X>,
        ::PMacc::cursor::BufferCursor<float_X, DIM1> >::type LinInterpCursor;

    typedef float_X type;

    LinInterpCursor linInterpCursor;

    /** constructor
     *
     * @param linInterpCursor lookup table of the first or the second
     * synchrotron function.
     */
    HDINLINE MapToLookupTable(LinInterpCursor linInterpCursor)
        : linInterpCursor(linInterpCursor) {}

    /** Returns F_1(x) or F_2(x)

     * @param x position of the synchrotron function to be evaluated
     */
    HDINLINE float_X operator()(const float_X x) const;
};

typedef ::PMacc::cursor::Cursor<
    MapToLookupTable,
    ::PMacc::cursor::PlusNavigator,
    float_X> SyncFuncCursor;

} // namespace detail


/** Lookup table for synchrotron functions.
 *
 * Provides cursors for the first and the second synchrotron function
 */
class SynchrotronFunctions
{
public:
    typedef detail::SyncFuncCursor SyncFuncCursor;
private:

    typedef boost::shared_ptr<PMacc::container::DeviceBuffer<float_X, DIM1> > MyBuf;
    MyBuf dBuf_SyncFuncs[2]; // two synchrotron functions

    struct BesselK
    {
        template<typename T_State, typename T_Time>
        void operator()(const T_State &x, T_State &dxdt, T_Time t) const
        {
            dxdt[0] = boost::math::tr1::cyl_bessel_k(5.0/3.0, t);
        }
    };

    /** First synchrotron function
     */
    float_64 F_1(const float_64 x) const;
    /** Second synchrotron function
     */
    float_64 F_2(const float_64 x) const;

public:
    enum Select
    {
        first=0, second=1
    };

    void init();
    /** Return a cursor representing a synchrotron function
     *
     * @param syncFunction first or second synchrotron function
     * @see: SynchrotronFunctions::Select
     */
    SyncFuncCursor getCursor(Select syncFunction) const;

}; // class SynchrotronFunctions

} // namespace synchrotronPhotons
} // namespace particles
} // namespace picongpu
