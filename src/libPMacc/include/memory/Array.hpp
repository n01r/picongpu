/**
 * Copyright 2016-2017 Rene Widera
 *
 * This file is part of libPMacc.
 *
 * libPMacc is free software: you can redistribute it and/or modify
 * it under the terms of either the GNU General Public License or
 * the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libPMacc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License and the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * and the GNU Lesser General Public License along with libPMacc.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "pmacc_types.hpp"


namespace PMacc
{
namespace memory
{
    /** static sized array
     *
     * mimic the most parts of the `std::array`
     */
    template<
        typename T_Type,
        size_t T_size
    >
    struct Array
    {
        using value_type = T_Type;
        using size_type = size_t;
        using reference = value_type &;
        using const_reference = value_type const &;
        using pointer = value_type *;
        using const_pointer = value_type const *;

        /** get number of elements */
        HDINLINE
        constexpr size_type size( ) const
        {
            return T_size;
        }

        /** get maximum number of elements */
        HDINLINE
        constexpr size_type max_size( ) const
        {
            return T_size;
        }

        /** get the direct access to the internal data
         *
         * @{
         */
        HDINLINE
        pointer data( )
        {
            return m_data;
        }

        HDINLINE
        const_pointer data( ) const
        {
            return m_data;
        }
        /** @} */

        /** default constructor
         *
         * the default constructor of each member is called
         */
        HDINLINE Array() = default;

        /** constructor
         *
         * initialize each member with the given value
         *
         * @param value element assigned to each member
         */
        HDINLINE Array( T_Type const & value )
        {
            for( size_type i = 0; i < size(); ++i )
                m_data[ i ] = value;
        }

        /** get N-th value
         *
         * @tparam T_Idx any type which can be implicit casted to an integral type
         * @param idx index within the array
         *
         * @{
         */
        template< typename T_Idx >
        HDINLINE
        const_reference
        operator[]( T_Idx const idx ) const
        {
            return m_data[ idx ];
        }

        template< typename T_Idx >
        HDINLINE
        reference
        operator[]( T_Idx const idx )
        {
            return m_data[ idx ];
        }
        /** @} */

    private:
        /** data storage */
        value_type m_data[ T_size ];
    };

} // namespace memory
} // namespace PMacc
