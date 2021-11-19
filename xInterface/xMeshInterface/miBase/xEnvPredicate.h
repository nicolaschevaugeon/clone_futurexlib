#ifndef _XENVPREDICATE_H_
#define _XENVPREDICATE_H_

#include "xEntity.h"


namespace xinterface{

  namespace xmeshinterface{


    class xEnvPredicate{
    public:
    xEnvPredicate(const xEnvPredicate &in) : iClassId(in.iClassId) , iClassLevel(in.iClassLevel) {};
    xEnvPredicate(int id ,int dim) : iClassId(id), iClassLevel(dim) {};
    private: 
      int iClassId, iClassLevel;
    public:
 
      inline  bool operator()(const xEntity &e){
     
       
	return ( !e.isAdjacencyCreated( e.getLevel() ) &&
		 e.getClassification().tag() == iClassId &&
		 e.getClassification().dim() == iClassLevel
		 );
           
      };
    };


  } // namepsace xmeshinterface
} // namespace xinterface

#endif
