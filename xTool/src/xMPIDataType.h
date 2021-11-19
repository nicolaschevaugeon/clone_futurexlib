/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _XMPIDATATYPE__H
#define _XMPIDATATYPE__H
#include <complex>
#include <iostream>
#include "mpi.h"
namespace xtool
{
////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// xMPIDataType function //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
 * This template function return the MPI data type corresponding to the template parameter type.
 * It works by specialization.
 */
template <typename T>
MPI_Datatype xMPIDataType()
{
   std::cerr << "Your type is not specialized !" << std::endl;
   throw - 5986341;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// xMPIDataType function specialization ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
 * In this location you will found only real generic types. Please do not add to much here (i.e. do
 * not add dependancy). The idea is that you add your own specialization (still in xtool namespace)
 * where your type is declared.
 */
template <>
MPI_Datatype xMPIDataType<int>();

template <>
MPI_Datatype xMPIDataType<size_t>();

template <>
MPI_Datatype xMPIDataType<float>();

template <>
MPI_Datatype xMPIDataType<double>();

template <>
MPI_Datatype xMPIDataType<std::complex<double>>();
}

#endif
