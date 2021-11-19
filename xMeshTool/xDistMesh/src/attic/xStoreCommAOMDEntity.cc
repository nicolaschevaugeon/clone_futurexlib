/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#include "xStoreCommAOMDEntity.h"
// std
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <ctime>

//AOMD
#include "mAOMD.h"
#include "mEntity.h"

using namespace std;
using namespace AOMD;


namespace distmesh
{
//===== xStoreCommEntity ===================================================================================================
xStoreCommEntity < AOMD::mEntity >::xStoreCommEntity(void) : proc_id(0), proc_id_is_assigned(false)
{
    stringstream stringtag;
    srand(time(NULL));
    stringtag << "XStoreCommEntitytagbnd" <<this<<rand();
    tag_bnd = AOMD::AOMD_Util::Instance()->lookupMeshDataId(stringtag.str().c_str());
}
xStoreCommEntity < AOMD::mEntity >::xStoreCommEntity(const typename xStoreCommEntity < AOMD::mEntity >::xStoreCommEntity & other) : proc_id(other.proc_id), proc_id_is_assigned(other.proc_id_is_assigned)
{
    stringstream stringtag;
    srand(time(NULL));
    stringtag << "XStoreCommEntitytagbnd" <<this<<rand();
    tag_bnd = AOMD::AOMD_Util::Instance()->lookupMeshDataId(stringtag.str().c_str());
    for (auto it = other.attached.begin(),ite = other.attached.end(); it != ite; ++it) 
    {
        AOMD::mEntity *e=*it;
        xCommEntity *ce = other.getAttachedCommEntity(e);
        attachCommEntity(e, *ce);
        attached.insert(e);
    }
}
xStoreCommEntity < AOMD::mEntity >::~xStoreCommEntity(void)
{
    for (auto it = attached.begin(),ite = attached.end(); it != ite; ++it) deleteCommEntity(*it);
    AOMD::AOMD_Util::Instance()->deleteMeshDataId (tag_bnd);
}
void xStoreCommEntity < AOMD::mEntity >::set(AOMD::mEntity *ladress, size_t rproc, AOMD::mEntity * radress)
{

    xCommEntity *ce = getAttachedCommEntity(ladress);
    void *ra = static_cast < void * >( radress );
    if (ce)
    {
        ce->add(xCommConector(rproc,ra));
    }
    else
    {
        xCommEntity c;
        c.add ( xCommConector(rproc,ra));
        attachCommEntity(ladress, c);
        attached.insert(ladress);
    }
    return;
}
xCommEntity * xStoreCommEntity < AOMD::mEntity >::get(AOMD::mEntity* e)
{
    return getAttachedCommEntity(e);
}
void xStoreCommEntity < AOMD::mEntity >::clear(AOMD::mEntity* e)
{
    if ( deleteCommEntity(e) )
    {
        attached.erase(e);
    }
}
void xStoreCommEntity < AOMD::mEntity >::clear(AOMD::mEntity* e, size_t rproc)
{
    xCommEntity *ce = getAttachedCommEntity(e);
    if (ce)
    {
        ce->remove(rproc);
        if (!ce->size()) clear(e);
    }
}
bool xStoreCommEntity < AOMD::mEntity >::deleteCommEntity(AOMD::mEntity* e)
{
    xAttachableCommEntity *av = (xAttachableCommEntity *) ( e->getData(tag_bnd));
    if (!av) return false;
    e->deleteData(tag_bnd);
    return true;
}
xCommEntity * xStoreCommEntity < AOMD::mEntity >::getAttachedCommEntity(AOMD::mEntity* e) const
{
    xAttachableCommEntity *av = (xAttachableCommEntity *) ( e->getData(tag_bnd));
    if (!av) return 0;
    return &( av->ace );
}
void xStoreCommEntity < AOMD::mEntity >::attachCommEntity(AOMD::mEntity* e, xCommEntity & v)
{
    xAttachableCommEntity*av = (xAttachableCommEntity *) ( e->getData(tag_bnd));
    if (!av)
    {
        av = new xAttachableCommEntity;
        e->attachData(tag_bnd,av);
    }
    av->ace = v;
}
void xStoreCommEntity < AOMD::mEntity >::assignProcId(size_t proc_id_)
{
    proc_id = proc_id_;
    proc_id_is_assigned = true;
}
bool xStoreCommEntity < AOMD::mEntity >::isOwner(AOMD::mEntity *e)
{
    xCommEntity* ce = getAttachedCommEntity(e);
    if (ce)
    {
        if (proc_id_is_assigned)
        {
            return ce->firstRemotIsGreaterThen(proc_id);
        }
        else
        {
            cout<<"You must assign a proc id before using isOwner method !"<<endl;
            throw -1;
        }
    }
    else
        return true;
}
//===== xAttachableCommEntity ==============================================================================================
xStoreCommEntity < AOMD::mEntity >::xAttachableCommEntity::~xAttachableCommEntity() {}

}     // end namspace
