#ifndef _MI_XCLASSIFICATION_H_
#define _MI_XCLASSIFICATION_H_


namespace xinterface{

  namespace xmeshinterface{

    class xClassification 
    {
    private:
      void* Geom; // = pointer to a geometric entity : equivalent to AOMD::pGEntity ou AOMD::GEntity*
      int Tag;
      int Dim;
    public:
    xClassification( void* _geom , int _tag, int _dim): Geom(_geom), Tag(_tag), Dim(_dim) {};
      ~xClassification(){};
      void* geom() const { return Geom; };
      int tag()    const { return Tag ; };
      int dim()    const { return Dim ; };
    };


  } // namepsace xmeshinterface
} // namespace xinterface

#endif
