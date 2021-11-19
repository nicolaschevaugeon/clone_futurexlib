/* 
   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of Trellis written and maintained by the 
   Scientific Computation Research Center (SCOREC) at Rensselaer Polytechnic
   Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA
*/
// mVector.cpp: implementation of the mVector class.
//
//////////////////////////////////////////////////////////////////////
#include "mVector.h"
#include "mTensor2.h"
#include "mPoint.h"

#include <cmath>

namespace Trellis_Util {

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

        mVector& mVector::operator *= ( const mTensor2 &other )
        {
            //      RIGHT multiply of a matrix to a vector: mTensor2 * mVector
            mVector m(*this);
            //flatten loops...
            pos[0] = other(0,0) * m.pos[0] + other(0,1) * m.pos[1] + other(0,2) * m.pos[2];
            pos[1] = other(1,0) * m.pos[0] + other(1,1) * m.pos[1] + other(1,2) * m.pos[2];
            pos[2] = other(2,0) * m.pos[0] + other(2,1) * m.pos[1] + other(2,2) * m.pos[2];

            return *this;
        }

        mVector mVector::operator * ( const mTensor2& M ) const
        {
        //      added by T. Bui 2/03
        //      LEFT multiply of a vector to a matrix: mVector * mTensor2
                return !M * (*this);
        }

} // end of namespace
