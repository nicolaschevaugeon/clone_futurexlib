/**************************************************************************** 

   Copyright (C) 2004
   Rensselaer Polytechnic Institute

   This file is part of the Algorithm-Oriented Mesh Database (AOMD) written 
   and maintained by the Scientific Computation Research Center (SCOREC) at 
   Rensselaer Polytechnic Intitute, Troy, NY, USA.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the Rensselaer SCOREC Public License.

   This program is distributed in the hope that it will be useful, 
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   license text for more details.
   
   You should have received a copy of the Rensselaer SCOREC Public License
   along with this program; if not, write to Rensselaer Polytechnic Institure,
   110 8th Street, SCOREC, Troy, NY  12180, USA

*****************************************************************************/
#include "mException.h"

#include <cstring>
#include <cstdio>

namespace AOMD {

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

mException::mException(const int i, const char *f, const char *text)
: theLine(i)
{
   strcpy(theFile,f);
   printf("an exception has been thrown at line %d in file %s\n",i,f);
   printf("TEXT : %s\n",text);
}

mException::~mException()
= default;

} // end of namespace
