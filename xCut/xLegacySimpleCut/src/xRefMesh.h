/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms 
   and conditions.
*/

#ifndef ___XREFMESH__H
#define ___XREFMESH__H

// All implemented things under this macro are not necesary for a On-the-fly cut during integrationRule
// They are implemented only for comptability whit old xPhysSurf
// They will be more simple to remove when a On-the-fly or equivalent implementation occures (undefine the macro or remove related portion of code)
//
// Default is to use sub definition so we echo a warning if not set
// For Now XREFMESH_WITH_SUB must be set, since otherwise XFEM library won't compil ... It needs to be tested if we wan't this option of not creating sub entities to work. Ideally it should be possible to set it dynamically. 
// Uncoment if you don't want this to be set by -DXREFMESH_WITH_SUB argument
#define XREFMESH_WITH_SUB 1
//
// Default is to use sub definition so we echo a warning if not set
#ifndef XREFMESH_WITH_SUB
#warning "In Xrefmesh compilation XREFMESH_WITH_SUB is not set"
#endif

#include <string>
#include <vector>
#include <map>

namespace xcut
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// xPoint struct //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  struct xPoint { 
  public: 
    double u;
    double v;
    double w;
  }; 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// End xPoint struct //////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// xRefMesh class /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class xRefMesh { 
  public: 
    xRefMesh(){};
    virtual ~xRefMesh()= default;

    // types
    typedef int idx_t;
    // container type
    typedef std::map<idx_t,xPoint> point_t;
    typedef std::map<idx_t,std::vector<idx_t> > elem_t;
    typedef std::vector<idx_t>  elemdef_t;
    // iterator
    typedef std::map<idx_t,xPoint>::iterator point_it;
    typedef std::map<idx_t,std::vector<idx_t> >::iterator elem_it;

    // public methodes ///////////////
    //
    // adding
    void addPoint(idx_t , xPoint & );
    void addElem(idx_t ,  elemdef_t & );
#ifdef XREFMESH_WITH_SUB
    void addEdge(idx_t ,  elemdef_t & );
    void addFace(idx_t ,  elemdef_t & );
#endif

    // getting
    xPoint & getPoint(idx_t);
    elemdef_t & getElem(idx_t);
    idx_t getLimite() {  return(limite); };
    idx_t getIsozero() {  return(iso_zero); };

    // setting
    void setLimite(idx_t id) { limite = id; };
    void setIsozero(idx_t id) { iso_zero = id; };
#ifdef XREFMESH_WITH_SUB
    void setEdgeLimite(idx_t id) { edge_limite = id; };
    void setFaceLimite(idx_t id) { face_limite = id; };
    void setIsozeroFace(idx_t id) { iso_zero_face = id; };
#endif

    // iterator
    elem_it elemBegin();
    elem_it elemEnd();
    elem_it elemLim();
    elem_it elemIsozero();
    point_it pointBegin();
    point_it pointEnd();
#ifdef XREFMESH_WITH_SUB
    elem_it edgeBegin();
    elem_it edgeEnd();
    elem_it edgeLim();
    elem_it faceBegin();
    elem_it faceEnd();
    elem_it faceLim();
    elem_it faceIsozero();
#endif

    // TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO
    // many over methode are needed to make this class realy user friendly
    // const and const iterator needed
    // comment needed to
    // TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO --- TODO
    // public members ////////////////
       
  private:
    // private members ////////////////
    point_t points;
    elem_t elems;
    idx_t limite;
    idx_t iso_zero;
#ifdef XREFMESH_WITH_SUB
    elem_t edges;
    idx_t edge_limite;
    elem_t faces;
    idx_t face_limite;
    idx_t iso_zero_face;
#endif
    // private methodes ///////////////


  };
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// End xRefMesh class /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// xRefMeshException class ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// interface derived class of standart exception for xsurf
  class xRefMeshException : public std::exception
  {
  public:
    xRefMeshException(std::string ,std::string ,int ,std::string ,std::string);
    ~xRefMeshException() throw() override;
    const char * what() const throw() override;

  private:
    std::string msg;
  };
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////// End xRefMeshException  class ///////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

} // end of namespace

#endif
