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

typedef struct SC_MemBlockList {
  void * d_block;
  void * d_nextBlock;
} SC_MemBlockList;

typedef struct SC_Allocator {
  int d_allocSize; /* size in bytes of memory being allocated (size of chunk) */
  int d_numInBlock; /* number of chunks to allocate at once */
  SC_MemBlockList *d_firstBlock;
  SC_MemBlockList *d_lastBlock;
  void * d_nextFree; /* next free location, null if none */
  int d_numAlloc;
  int d_numFree;
} SC_Allocator;

SC_Allocator * SC_newAllocator(int allocSize, int blockSize);
void SC_deleteAllocator(SC_Allocator *);
void * SC_allocate(SC_Allocator *);
void SC_deallocate(SC_Allocator *b, void *mem);
