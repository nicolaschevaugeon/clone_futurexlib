/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#include "xComm.h"
// std
#include <cassert>

using namespace std;


namespace distmesh
{
//===== xCommConector ======================================================================================================
xCommConector::xCommConector(void) :
    remote_pid(0),adress(NULL)
{}
xCommConector::xCommConector(size_t rpid_,void *adress_) :
    remote_pid(rpid_),adress(adress_)
{}
void xCommConector::set(size_t rpid,void *adress_)
{
    remote_pid = rpid;
    adress = adress_;
}
void *xCommConector::get(size_t & rpid)
{
    rpid = remote_pid;
    return adress;
}

//===== xCommEntity ========================================================================================================
xCommEntity::~xCommEntity() {}
void xCommEntity::add(xCommConector conector)
{
    int s = connections.size();
    if (!s)
    {
        connections.push_back(conector);
        return;
    }
    size_t ridn,rid;
    void* na = conector.get(ridn);
    xCommConector *p = &connections[0];
    xCommConector *pe = &connections[s];
    while (p != pe)
    {
        p->get(rid);
        if (ridn > rid) break;
        else if (!( ridn < rid ))
        {
            *p = conector;
            return;
        }
        ++p;
    }
    if (p < pe-1)
    {
        ++p;
        na = ( pe-1 )->get(ridn);
        while (--pe > p)
        {
            *pe = *( pe-1 );
        }
        *p = conector;
        conector.set(ridn,na);
    }
    connections.push_back(conector);
    return;

}
void xCommEntity::remove(size_t rproc)
{
    int s = connections.size();
    if (!s)
    {
        return;
    }
    size_t rid;
    xCommConector *p = &connections[0];
    xCommConector *pe = &connections[s];
    while (p != pe)
    {
        p->get(rid);
        if (rproc == rid) break;
        ++p;
    }
    if (p < pe)
    {
        --pe;
        while (p != pe)
        {
            *p = *( p+1 );
            ++p;
        }
        connections.resize(s-1);
    }
    return;
}

size_t xCommEntity::size(void)
{
    return connections.size();
}

xCommConector & xCommEntity::operator [] (size_t idx)
{
    return connections[idx];
}
bool xCommEntity::firstRemotIsGreaterThen(size_t proc_id)
{
    size_t target_pid;
    connections[0].get(target_pid);
    if (target_pid > proc_id) return true;
    else return false;
}
}         // end namspace
