/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
 */
#ifndef XSPARSEVECTORPACKUNPACK_H
#define XSPARSEVECTORPACKUNPACK_H

// xlinalg
#include "xSparseVector.h"


// exlibris_tools
#include "xDataExchanger.h"
#include "xMPIDataType.h"

/// A function to pack
template <typename T>
void pack(const xlinalg::xSparseVector<T> &vect,xtool::xMpiInputBuffer &buff);

template <typename T>
void unPack(xlinalg::xSparseVector<T> &vect, const xtool::xMpiOutputBuffer &buff);

#include "xSparseVectorPackUnPack_imp.h"
#endif
