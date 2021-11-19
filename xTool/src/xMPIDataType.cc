/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
 */

#include "xMPIDataType.h"
#include <type_traits>

namespace xtool
{
////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// xMPIDataType function specialization implementation ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
template <>
MPI_Datatype xMPIDataType<int>()
{
   return MPI_INT;
}
template <>
MPI_Datatype xMPIDataType<size_t>()
{
   // from most to less probable
   // Anyway, hopefully, the compiler optimize this
   // part by removing all test as compared type are
   // all known at compilation time
   if (std::is_same<unsigned long, size_t>::value) return MPI_UNSIGNED_LONG;
   if (std::is_same<unsigned long long, size_t>::value) return MPI_UNSIGNED_LONG_LONG;
   if (std::is_same<unsigned int, size_t>::value) return MPI_UNSIGNED;
   if (std::is_same<unsigned short, size_t>::value) return MPI_UNSIGNED_SHORT;
   if (std::is_same<unsigned char, size_t>::value) return MPI_UNSIGNED_CHAR;
   throw - 3345;
}

template <>
MPI_Datatype xMPIDataType<float>()
{
   return MPI_REAL;
}

template <>
MPI_Datatype xMPIDataType<double>()
{
   return MPI_DOUBLE;
}

template <>
MPI_Datatype xMPIDataType<std::complex<double>>()
{
   return MPI_DOUBLE_COMPLEX;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// End xMPIDataType function specialization implementation ////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
}  // end of namespace
