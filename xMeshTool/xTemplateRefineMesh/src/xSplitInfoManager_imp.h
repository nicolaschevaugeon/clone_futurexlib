/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE & LICENSE files for terms
   and conditions.
 */

#ifndef XSPLITINFOMANAGER_IMP_H
#define XSPLITINFOMANAGER_IMP_H

namespace xmeshtool
{

//===== splitKeyManagerSendAndReceive ======================================================================================================
template < typename PM, typename T >
splitKeyManagerSendAndReceive < PM,T >::splitKeyManagerSendAndReceive( PM &partman_ ) : partman(partman_){};

template < typename PM, typename T >
auto splitKeyManagerSendAndReceive < PM,T >::localObjectKey( const T &o)->information_key_t
{
    return &o;
}

template < typename PM, typename T >
auto splitKeyManagerSendAndReceive < PM,T >::remoteObjectKey(const xtool::xRemoteObject < T > &ro, const T &lo )->information_key_t
{
    return ro.getObjectAddress();
}
template < typename PM, typename T >
xtool::xConstPartitionObject < T > splitKeyManagerSendAndReceive < PM,T >::getConstPartitionObject(const T &e)
{
    return partman.getConstPartitionObject(e);
}
//===== splitKeyManagerSendOnly ======================================================================================================
template < typename PM, typename T >
splitKeyManagerSendOnly < PM,T >::splitKeyManagerSendOnly( PM &partman_ ) : partman(partman_){};
template < typename PM, typename T >
auto splitKeyManagerSendOnly < PM,T >::localObjectKey( const T &o)->information_key_t
{
    return &o;
}
template < typename PM, typename T >
std::set < int > splitKeyManagerSendOnly < PM,T >::getMessageRanks( const information_key_t &lk)
{
    // lk is the pointer to local object
    T &lo = *( const_cast < T * >( lk ));
    // So we can used it to access to remote object
    xtool::xConstPartitionObject < T > po = partman.getConstPartitionObject(lo);

    // object accumulated are supposed to be bnd entity having hanging edge or face.
    // if they have no remote it's a problem
    assert (po.hasRemoteObject());

    // grabe remote proc id for object and store them in res
    std::set < int > res;
    for (const auto ro : po.getRemoteObjectsCollectionRange( ))
    {
        res.insert(ro.getProcId());
    }


    return res;
}
//===== splitInfoManagerOneLevelInformation ======================================================================================================
template < typename PM, typename T >
splitInfoManagerOneLevelInformation < PM,T >::splitInfoManagerOneLevelInformation(const PM &part_man_,std::set < T * > &to_be_added_) :
    partman( part_man_ )
    ,to_be_added(to_be_added_)
{}
template < typename PM, typename T >
auto splitInfoManagerOneLevelInformation < PM,T >::getInfo( information_key_t lk, int sendto )->information_t
{
    // lk is the pointer to local object
    // So we can used it to access to remote object
    xtool::xConstPartitionObject < T > po = partman.getConstPartitionObject(*lk);

    // object accumulated are supposed to be bnd entity having a one level rule problem
    // if they have no remote it's a problem
    assert (po.hasRemoteObject());

    // get remote object information
    auto ro = po.getRemoteObjectOn(sendto);

    // in getMessageRanks we did the job of giving remote proc id for this object
    // so there is no reason to have a null pointer
    assert(ro);

    return const_cast < information_t >( ro );
}
template < typename PM, typename T >
void splitInfoManagerOneLevelInformation < PM,T >::setInfo( const std::vector < information_t > &info, int receivedfrom )
{
    to_be_added.insert(info.cbegin(),info.cend());
}
//===== splitInfoManagerNewHangingTreatement ===============================================================================================================
template < typename PM, typename T >
splitInfoManagerNewHangingTreatement < PM,T >::splitInfoManagerNewHangingTreatement(const PM &part_man_,std::function < void(T &) > &treat_new_hanging_) :
    partman ( part_man_ )
    ,treat_new_hanging(treat_new_hanging_)
{}
template < typename PM, typename T >
auto splitInfoManagerNewHangingTreatement < PM,T >::getInfo( information_key_t lk, int sendto )->information_t
{
    // lk is the pointer to local object
    // So we can used it to access to remote object
    xtool::xConstPartitionObject < T > po = partman.getConstPartitionObject(*lk);

    // object accumulated are supposed to be new bnd entity having hanging edge or face. But these are not yet been synchronise with
    // remotes (i.e. update of partman for hanging and eventual spliting in remote proc)
    // if they have no remote it's a problem
    assert (po.hasRemoteObject());

    // get remote object information
    auto ro = po.getRemoteObjectOn(sendto);

    // in getMessageRanks we did the job of giving remote proc id for this object
    // so there is no reason to have a null pointer
    assert(ro);

    // this is to treate local object
    T &lo = *( const_cast < T * >( lk ));
    treat_new_hanging(lo);

    return const_cast < information_t >( ro );
}
template < typename PM, typename T >
void splitInfoManagerNewHangingTreatement < PM,T >::setInfo( const std::vector < information_t > &info, int receivedfrom )
{
    // this is to treat remote object
    for (const auto inf : info) treat_new_hanging(*inf);
}

//===== splitInfoManagerUpdatePMFE   ===============================================================================================================
template < typename PM, typename T, typename DM, typename DM1, typename DM2 >
splitInfoManagerUpdatePMFE < PM,T,DM,DM1,DM2 >::splitInfoManagerUpdatePMFE( PM &part_man_,DM & hd, DM & hu, DM1 &st1, DM2 &st2) :
    partman(part_man_)
    ,hanging_down(hd)
    ,hanging_up(hu)
    ,sub_t1(st1)
    ,sub_t2(st2)
{}
template < typename PM, typename T, typename DM, typename DM1, typename DM2 >
void splitInfoManagerUpdatePMFE < PM,T,DM,DM1,DM2 >::getInfo(information_key_t key, xtool::xMpiInputBuffer & buff, int sendto)
{
    T * const *ppe = hanging_down.getData(*key);
    if (ppe)
    {
        T *pe = *ppe;
        // for now we have only face or node respectively for face and edge. No switch, just a test

        // face with hanging
        if (pe->getLevel() > 0)
        {
            assert(pe->size(1) == 3);
            buff.pack(&pe,1,MPI_AINT);
            bool found;
            for (int i = 0,n = 3; i < n; ++i)
            {
                T *edge = pe->get(1,i);
                buff.pack(&edge,1,MPI_AINT);
                assert(edge->size(2));
                int j,m;
                found = false;
                for (j = 0,m = edge->size(2); j < m; ++j)
                {
                    T *face = edge->get(2,j);
                    if (face != pe)
                    {
                        T **puface = hanging_up.getData(*face);
                        if (puface && *puface == key)
                        {
                            buff.pack(&face,1,MPI_AINT);
                            found = true;
                            break;
                        }
                    }
                }
                assert(found);
            }
            found |= true; //for -Wunused-but-set-variable
        }
        // edge with hanging
        else
        {
            T & mid_node = *pe;
            assert(mid_node.size(1));
            buff.pack(&pe,1,MPI_AINT);
            T *v1 = key->get(0,0);
            std::array < T *, 2 >  info;
            int m = 0;
            for (int i = 0,n = mid_node.size(1); i < n && m < 2; ++i)
            {
                T *edge = mid_node.get(1,i);
                T **puedge = hanging_up.getData(*edge);
                if (puedge && *puedge == key)
                {
                    if (edge->get ( 0,0) == v1)
                    {
                        info[0] = edge;
                        ++m;
                    }
                    else
                    {
                        info[1] = edge;
                        ++m;
                    }
                }
            }
            assert( m == 2);
            buff.pack(info.data(),2,MPI_AINT);
        }
    }
    else
    {
        // use of auto to mask real pointer to array stuff
        auto related_face_t1 = sub_t1.getData(*key);
        // if stored in sub_t1
        if (related_face_t1)
            buff.pack(related_face_t1->data(),3,MPI_AINT);
        else
        {
            // use of auto to mask real pointer to array stuff
            auto related_face_t2 = sub_t2.getData(*key);
            // if stored in sub_t2
            if (related_face_t2)
                buff.pack(related_face_t2->data(),5,MPI_AINT);
            //it should be on of those type. if not it is an error
            else
            {
                throw -12345;
            }
        }
    }

    return;

}
template < typename PM, typename T, typename DM, typename DM1, typename DM2 >
void splitInfoManagerUpdatePMFE < PM,T,DM,DM1,DM2 >::setInfo(information_key_t key, const xtool::xMpiOutputBuffer & buff, int receivedfrom)
{
    T * const *ppe = hanging_down.getData(*key);
    if (ppe)
    {
        T *pe = *ppe;
        T *tmp;
        // for now we have only face or node respectively for face and edge. No switch, just a test

        // face with hanging
        if (pe->getLevel() > 0)
        {
            xtool::xPartitionObject < T, PM > po = partman.getPartitionObject(*pe);
            buff.unPack(&tmp,1,MPI_AINT);
            po.insert(receivedfrom,tmp);
            assert(pe->size(1) == 3);
            for (int i = 0,n = 3; i < n; ++i)
            {
                T &edge = *( pe->get(1,i));
                xtool::xPartitionObject < T, PM > poi = partman.getPartitionObject(edge);
                buff.unPack(&tmp,1,MPI_AINT);
                poi.insert(receivedfrom,tmp);
                assert(edge.size(2));
                int j,m;
                for (j = 0,m = edge.size(2); j < m; ++j)
                {
                    T *face = edge.get(2,j);
                    if (face != pe)
                    {
                        T **puface = hanging_up.getData(*face);
                        if (puface && *puface == key)
                        {
                            xtool::xPartitionObject < T, PM > poij = partman.getPartitionObject(*face);
                            buff.unPack(&tmp,1,MPI_AINT);
                            poij.insert(receivedfrom,tmp);
                            break;
                        }
                    }
                }
                assert(j < m);
            }
        }
        // edge with hanging
        else
        {
            T & mid_node = *pe;
            assert(mid_node.size(1));
            std::array < T *, 3 >  info;
            buff.unPack(info.data(),3,MPI_AINT);
            xtool::xPartitionObject < T, PM > po = partman.getPartitionObject(mid_node);
            po.insert(receivedfrom,info[0]);

            T *v1 = key->get(0,0);
            int m = 0;
            for (int i = 0,n = mid_node.size(1); i < n && m < 2; ++i)
            {
                T &edge = *( mid_node.get(1,i));
                T **puedge = hanging_up.getData(edge);
                if (puedge && *puedge == key)
                {
                    xtool::xPartitionObject < T, PM > poi = partman.getPartitionObject(edge);
                    if (edge.get(0,0) == v1)
                    {
                        poi.insert(receivedfrom,info[1]);
                        ++m;
                    }
                    else
                    {
                        poi.insert(receivedfrom,info[2]);
                        ++m;
                    }
                }
            }
            assert( m == 2);
        }
    }
    else
    {
        // use of auto to mask real pointer to array stuff
        auto related_face_t1 = sub_t1.getData(*key);
        // if stored in sub_t1
        if (related_face_t1)
        {
            T *tmp;
            for (auto enti :*related_face_t1)
            {
                buff.unPack(&tmp,1,MPI_AINT);
                xtool::xPartitionObject < T, PM > po = partman.getPartitionObject(*enti);
                po.insert(receivedfrom,tmp);
            }
        }
        else
        {
            // use of auto to mask real pointer to array stuff
            auto related_face_t2 = sub_t2.getData(*key);
            // if stored in sub_t2
            if (related_face_t2)
            {
                T *tmp;
                for (auto enti :*related_face_t2)
                {
                    buff.unPack(&tmp,1,MPI_AINT);
                    xtool::xPartitionObject < T, PM > po = partman.getPartitionObject(*enti);
                    po.insert(receivedfrom,tmp);
                }
            }
            //it should be on of those type. if not it is an error
            else
            {
                throw -12874;
            }
        }
    }
    return;
}
template < typename PM, typename T, typename DM, typename DM1, typename DM2 >
size_t splitInfoManagerUpdatePMFE < PM,T,DM,DM1,DM2 >::getApproxDataSize()
{
    return 3*sizeof( MPI_AINT );
}
//===== splitInfoManagerrDanglingCount ===============================================================================================================
template <  typename T, template < typename >  class DM >
void splitInfoManagerDanglingCount < T, DM >::setNbConnect(const information_key_t & k, const information_t &inf)
{
    nb_component.setData(*const_cast < T * >( k )) = inf;
}
template <  typename T, template < typename >  class DM >
void splitInfoManagerDanglingCount < T,DM >::getNbConnect(const information_key_t & k, information_t &inf)
{
    const information_t *pinf = nb_component.getData(*k);
    assert(pinf);
    inf = *pinf;
}
template <  typename T, template < typename >  class DM >
void splitInfoManagerDanglingCount < T,DM >::clearNbConnect(const information_key_t & k)
{
    nb_component.deleteData(*const_cast < T * >( k ));
}
template <  typename T, template < typename >  class DM >
auto splitInfoManagerDanglingCount < T,DM >::getInfo( information_key_t lk, int sendto )->information_t
{
    const information_t *pinf = nb_component.getData(*lk);
    assert(pinf);
    return *pinf;
}
template <  typename T, template < typename >  class DM >
void splitInfoManagerDanglingCount < T,DM >::setInfo(information_key_t lk, const information_t &info, int receivedfrom)
{
    const information_t *pinf = nb_component.getData(*lk);
    assert(pinf);
    nb_component.setData(*const_cast < T * >( lk )) = std::max(*pinf,info);
}

} //end namespace
#endif
