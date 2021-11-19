#include "SC_mem.h"
#include <cstdlib>
#include <cstdio>

/* note: allocSize must be a multiple of sizeof(void*)
   */

void *SC_allocNewBlock(int allocSize, int numAlloc)
     /* allocate the new block of memory and set up pointers for each
	allocatable piece in the new block to point to the next 
	available block */
{
  int i,allocInc;
  void **block, **bp;
  allocInc = allocSize/sizeof(void*);
#ifdef DEBUG
#ifndef SGI
  allocInc++; /* to put in guard block */
  allocSize+=sizeof(void*);
#endif
#endif
  block = (void**)malloc(allocSize*numAlloc);
  if(block){
    bp = block;
    for(i=0; i < numAlloc-1; i++){
      *(bp) = bp+allocInc; /* store pointer to next block */
#ifdef DEBUG
#ifndef SGI
      *(long*)(bp+allocInc-1) = allocSize-sizeof(void*);
#endif
#endif
      bp += allocInc;
    }
    *(bp) = nullptr; /* last one*/
#ifdef DEBUG
#ifndef SGI
    *(long*)(bp+allocInc-1) = allocSize-sizeof(void*);
#endif
#endif

  }
  return block;
}

SC_Allocator * SC_newAllocator(int allocSize, int numInBlock)
{
  SC_Allocator *block;
  block = (SC_Allocator*)malloc(sizeof(SC_Allocator));
  block->d_allocSize = allocSize;
  block->d_numInBlock = numInBlock;
  block->d_firstBlock = (SC_MemBlockList*)malloc(sizeof(SC_MemBlockList));
  block->d_firstBlock->d_nextBlock = nullptr;
  block->d_firstBlock->d_block = SC_allocNewBlock(allocSize,numInBlock);
  block->d_nextFree = block->d_firstBlock->d_block;

  block->d_lastBlock = block->d_firstBlock;
  block->d_numAlloc = numInBlock;
  block->d_numFree = numInBlock;
  return block;
}

void SC_deleteAllocator(SC_Allocator *b)
{
  SC_MemBlockList *bl, *bln;
  bl = b->d_firstBlock;
  while(bl){
    free(bl->d_block);
    bln = (SC_MemBlockList*)bl->d_nextBlock;
    free(bl);
    bl = bln;
  }
}

void * SC_allocate(SC_Allocator *b)
{
  void *ret;
  SC_MemBlockList *bl;
  ret = b->d_nextFree;
  if(!ret){
    bl = (SC_MemBlockList*)malloc(sizeof(SC_MemBlockList));
    bl->d_block = SC_allocNewBlock(b->d_allocSize,b->d_numInBlock);
    bl->d_nextBlock = nullptr;
    b->d_lastBlock->d_nextBlock = bl;
    b->d_lastBlock = (SC_MemBlockList*)b->d_lastBlock->d_nextBlock;
    b->d_nextFree = bl->d_block;
    ret = b->d_nextFree;
    b->d_numAlloc += b->d_numInBlock;
    b->d_numFree += b->d_numInBlock;
  } 
  b->d_nextFree = *((void**)(b->d_nextFree));
  b->d_numFree--;
  if(b->d_numFree < 0)
    fprintf(stderr, "allocator error 1\n");
#ifdef DEBUG
#ifndef SGI
  if(*( ((long*)ret)+b->d_allocSize/sizeof(void*)) != b->d_allocSize)
    fprintf(stderr,"allocate error 3\n");
#endif
#endif
  return ret;
}

void SC_deallocate(SC_Allocator *b, void *mem)
{
#ifdef DEBUG
#ifndef SGI
  if(*( ((long*)mem)+b->d_allocSize/sizeof(void*)) != b->d_allocSize)
    fprintf(stderr,"allocate error 4\n");
#endif
#endif
  *((void**)(mem)) = (void*)b->d_nextFree;
  b->d_nextFree = mem;
  b->d_numFree++;
}

void SC_allocatorState(SC_Allocator *b)
{
  fprintf(stderr,"allocSize = %d\n",b->d_allocSize);
  fprintf(stderr,"%d allocated\n",b->d_numAlloc);
  fprintf(stderr,"%d free\n",b->d_numFree);
}
