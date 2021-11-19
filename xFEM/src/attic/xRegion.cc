/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <string>
#include <iostream>
#include <cstdio>
#include "xRegion.h"
#include "xSubMesh.h"
#include "xMesh.h"

#ifdef PARALLEL 
#include "AOMD_OwnerManager.h"
#endif



namespace xfem
{
 
xRegion::xRegion() : mesh(0), sub(0), dimension(-1) {}
xRegion::xRegion(const xMesh* m) : mesh(const_cast<xMesh*>(m)), sub(0), dimension(m->dim()) {}
  xRegion::xRegion(const xMesh* m, const std::string& s) : mesh(const_cast<xMesh*>(m)), sub( &m->getSubMesh(s)), dimension(m->dim()) {}

  xRegion::xRegion(const xSubMesh * ms) : mesh(const_cast<xMesh*>( &ms->getMesh())), sub( ms), dimension(mesh->dim()) {}

xMesh* xRegion::getMesh() const { return mesh; }
std::string xRegion::getSubName() const { 
  if (sub) return sub->getName();
  else return std::string("");
}

void xRegion::setMesh(xMesh *m){
  mesh = const_cast<xMesh*>(m);
  //Faire fct le sub si besoin:
  //sub(0)
  dimension = m->dim();
}


xIter xRegion::begin() const {
  if (sub) return sub->begin(dimension);
  else return mesh->begin(dimension);
}
xIter xRegion::end() const {
  if (sub) return sub->end(dimension);
  else return mesh->end(dimension);
}

int  xRegion::size() const {
  if (sub) return sub->size(dimension);
  else return mesh->size(dimension);
}

xIter xRegion::begin(int what) const {
  if (sub) return sub->begin(what);
  else return mesh->begin(what);

}
xIter xRegion::end(int what) const {
  if (sub) return sub->end(what);
  else return mesh->end(what);
}

int  xRegion::size(int what) const {
  if (sub) return sub->size(what);
  else return mesh->size(what);
}

bool xRegion::IsInRegion(AOMD::mEntity* e) const {
  if (sub)return (sub->find(e))?true:false;
  else return (mesh->find(e))?true:false;
}
  
AOMD::mEntity* xRegion::find(AOMD::mEntity*e) const { 
  if (sub) return sub->find(e);
  else  return mesh->find(e);
}
  
#ifdef PARALLEL 
// get Parallel partition boundary information container
// it give the container from submesh or mesh if submesh don't exist
// nota : getAllSharedInfo is a addon of "AOMD_OwnerManager.h". It have been commited 22/11/10. 
//        Be shure to use the correct version of AOMD
const xPartitionBoundaryInfo & xRegion::getPartitionBoundaryInfo(void) const
{
   if (sub) return sub->getPartitionBoundaryInfo();
   else return( (mesh->theOwnerManager)->getAllSharedInfo());
}
#endif
// 

//  bool xRegion::operator==(const xRegion &other){
//    if (mesh != other.mesh)           return false;
//    if 
//    switch(iter_type) {
//    case BASIC_ITERATOR: return (sub == other.sub); break;
//    case CLASS_ITERATOR: return (entity_id == other.entity_id 
//  			       && entity_dim == other.entity_dim); break; 
//    case LEAVES_ITERATOR: assert(0); break;
//    }
//  }



xClassRegion::xClassRegion(xMesh* m, int e_id, int e_dim) : 
   mesh(m), entity_id(e_id), entity_dim(e_dim), dimension(e_dim) {}

xClassIter xClassRegion::begin() {
  return mesh->begin(entity_dim, entity_id, entity_dim);
}
xClassIter xClassRegion::end() {
  return mesh->end(entity_dim, entity_id, entity_dim);
}

xClassIter xClassRegion::begin(int what) {
  return mesh->begin(what, entity_id, entity_dim);
}
xClassIter xClassRegion::end(int what) {
  return mesh->end(what, entity_id, entity_dim);
}




} // end of namespace
