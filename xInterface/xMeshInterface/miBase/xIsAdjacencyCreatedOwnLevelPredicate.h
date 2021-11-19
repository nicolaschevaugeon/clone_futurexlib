#ifndef _XISADJACENCYCREATEDOWNLEVEL_PREDICATE_H_
#define _XISADJACENCYCREATEDOWNLEVEL_PREDICATE_H_

#include "xEntity.h"


namespace xinterface{

  namespace xmeshinterface{

    class xIsAdjacencyCreatedOwnLevelPredicate
    {
    public:
      xIsAdjacencyCreatedOwnLevelPredicate();
    public:
 
      inline  bool operator()(const xEntity &e)
      {
	return ( !e.isAdjacencyCreated( e.getLevel() ) );
      };
    };


  
  } // namepsace xmeshinterface
} // namespace xinterface
#endif
