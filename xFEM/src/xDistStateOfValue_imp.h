/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef XDISTSTATEOFVALUE_IMP_H
#define XDISTSTATEOFVALUE_IMP_H


namespace xfem
{
//===== finalizeFixedUniform =========================================================================================================
template < typename PM, typename VT >
bool finalizeFixedUniform(  const PM &part_man,  xValueManagerDist<VT> & dm)
{
    // create key manager
    xDistFixedKeyManagerUniform < PM, VT > key_info_manager(part_man, dm);

    // set keys for information to be exchanged
    xtool::xKeyContainerSendOrRecv < typename  xDistFixedKeyManagerUniform < PM, VT >::information_key_t, xHashValKey, xEqualValKey > infokey_container(part_man.getComm());
    infokey_container.accumulateKeys(dm.begin(),dm.end(),key_info_manager);

    // create info manager
    xDistFixedInfoManagerUniform < PM, VT > info_manager(part_man, dm);

    // exchange information (create fixed status if needed)
    exchangeInformation(infokey_container,info_manager);

    return ( info_manager.getNumberOfModification() != 0 );
}

//===== finalizeDof ===========================================================================================================
template < typename Condition, typename VT >
std::shared_ptr < xlinalg::xDistIndex > finalizeDof( xDistStateDofCreator < Condition, VT > &dof_creator)
{
    const typename xValueManagerDist<VT>::partman_t &part_man = dof_creator.getPartitionManager();
    auto &dm = dof_creator.getValueManager();
    const std::string sub = dof_creator.getSubset();
    MPI_Comm world = part_man.getComm();

    // create distributed index structure
    xlinalg::xDistIndex *dist_index = new xlinalg::xDistIndex(world);
    std::shared_ptr < xlinalg::xDistIndex > dist_index_shr_ptr (dist_index);

    int nbdof = dm.size(sub);
    const bool verbose = true;

    // create key(info key) and info manager
    xDistDofKeyAndInformationManager<VT> key_info_manager(dm,part_man);

    // set index of information to be exchanged
    xtool::xKeyContainerSendAndRecv < typename xDistDofKeyAndInformationManager<VT>::information_key_t > infokey_container(world);
    infokey_container.accumulateKeysOwnerScatter(dm.begin(),dm.end(),key_info_manager);

    if (verbose)
        std::cout<<"nbdof local "<<nbdof<< std::endl;

    // exchanging number of dof computed so far to create offset for renumbering
    int offset = 0;
    MPI_Scan ( &nbdof, &offset, 1, MPI_INT, MPI_SUM, world );
    int numdof_global = offset;
    offset -= nbdof;

    // renumbering local dof with offset and insert index
    if (offset)
    {
        for (auto it = dm.begin(sub),ite = dm.end(sub); it != ite; ++it)
        {
            xfem::xStateOfValueDof*state = static_cast < xfem::xStateOfValueDof * >(( *it )->getState());
            int n = state->Numdof + offset;
            dist_index->insertIndex(n);
            state->Numdof = n;
        }
    }
    // or just insert index
    else
    {
        for (auto it = dm.begin(sub),ite = dm.end(sub); it != ite; ++it)
        {
            xfem::xStateOfValueDof*state = static_cast < xfem::xStateOfValueDof * >(( *it )->getState());
            dist_index->insertIndex(state->Numdof);
        }
    }

    // state dof creator without distributed test
    xStateDofCreator < Condition, VT > cdof(dm, sub);

    // loop to create dumy dof before exchange
    for (auto it = dm.begin(),ite = dm.end(); it != ite; ++it)
    {
        cdof(it->first,dm.getValPtr(it->first));
    }

    // exchange common dof
    exchangeInformation(infokey_container,key_info_manager);

    // insert index for common dof
    for (auto it = dm.begin(sub)+nbdof,ite = dm.end(sub); it != ite; ++it)
    {
        xfem::xStateOfValueDof*state = static_cast < xfem::xStateOfValueDof * >(( *it )->getState());
        dist_index->insertIndex(state->Numdof);
    }

    // set remote duplicate index information
    for (auto it = dm.begin(),ite = dm.end(); it != ite; ++it)
    {
        auto po = key_info_manager.getConstPartitionObject(*it);
        if (po.hasRemoteObject())
        {
            xfem::xStateOfValue *state = it->second->getState();
            assert(state);
            if (state->state == xfem::xStateOfValue::DOF)
            {
                xfem::xStateOfValueDof *dof = static_cast < xfem::xStateOfValueDof * >( state );
                const int idx = dof->Numdof;
                for (const auto ro : po.getRemoteObjectsCollectionRange( )) dist_index->insertToFrom(idx,ro.getProcId());
            }
        }
    }

    // finalize index structure
    dist_index->finalize(nbdof,numdof_global,true);

    if (verbose)
        std::cout<<"nbdof global "<<dist_index->getGlobalIndexSize()<<std::endl;

    return dist_index_shr_ptr;
}
//===== finalizeDofUniform ========================================================================================================
template < typename PM, typename Condition, typename VT >
std::shared_ptr < xlinalg::xDistIndex > finalizeDofUniform( xDistStateDofCreatorUniform < PM, Condition, VT > &dof_creator)
{
    const PM &part_man = dof_creator.getPartitionManager();
    auto &dm = dof_creator.getValueManager();
    const std::string sub = dof_creator.getSubset();
    MPI_Comm world = part_man.getComm();

    // create distributed index structure
    xlinalg::xDistIndex *dist_index = new xlinalg::xDistIndex(world);
    std::shared_ptr < xlinalg::xDistIndex > dist_index_shr_ptr (dist_index);

    int nbdof = dm.size(sub);
    const bool verbose = true;

    // create key(info key) and info manager
    xDistDofKeyAndInformationManagerUniform < PM, VT > key_info_manager(part_man, dm);

    // set index of information to be exchanged
    xtool::xKeyContainerSendAndRecv < typename xDistDofKeyAndInformationManagerUniform < PM, VT >::information_key_t > infokey_container(world);
    infokey_container.accumulateKeysOwnerScatter(dm.begin(),dm.end(),key_info_manager);

    if (verbose)
        std::cout<<"nbdof local "<<nbdof<<std::endl;

    // exchanging number of dof computed so far to create offset for renumbering
    int offset = 0;
    MPI_Scan ( &nbdof, &offset, 1, MPI_INT, MPI_SUM, world );
    int numdof_global = offset;
    offset -= nbdof;

    // renumbering local dof with offset and insert index
    if (offset)
    {
        for (auto it = dm.begin(sub),ite = dm.end(sub); it != ite; ++it)
        {
            xfem::xStateOfValueDof*state = static_cast < xfem::xStateOfValueDof * >(( *it )->getState());
            int n = state->Numdof + offset;
            dist_index->insertIndex(n);
            state->Numdof = n;
        }
    }
    // or just insert index
    else
    {
        for (auto it = dm.begin(sub),ite = dm.end(sub); it != ite; ++it)
        {
            xfem::xStateOfValueDof*state = static_cast < xfem::xStateOfValueDof * >(( *it )->getState());
            dist_index->insertIndex(state->Numdof);
        }
    }

    // state dof creator without distributed test
    xStateDofCreator < Condition, VT > cdof(dm, sub);

    // loop to create dumy dof before exchange
    for (auto it = dm.begin(),ite = dm.end(); it != ite; ++it)
    {
        cdof(it->first,dm.getValPtr(it->first));
    }

    // exchange common dof
    exchangeInformation(infokey_container,key_info_manager);

    // insert index for common dof
    for (auto it = dm.begin(sub)+nbdof,ite = dm.end(sub); it != ite; ++it)
    {
        xfem::xStateOfValueDof*state = static_cast < xfem::xStateOfValueDof * >(( *it )->getState());
        dist_index->insertIndex(state->Numdof);
    }

    // set remote duplicate index information
    for (auto it = dm.begin(),ite = dm.end(); it != ite; ++it)
    {
        auto po = key_info_manager.getConstPartitionObject(*it);
        if (po.hasRemoteObject())
        {
            xfem::xStateOfValue *state = it->second->getState();
            assert(state);
            if (state->state == xfem::xStateOfValue::DOF)
            {
                xfem::xStateOfValueDof *dof = static_cast < xfem::xStateOfValueDof * >( state );
                const int idx = dof->Numdof;
                for (const auto ro : po.getRemoteObjectsCollectionRange( )) dist_index->insertToFrom(idx,ro.getProcId());
            }
        }
    }

    // finalize index structure
    dist_index->finalize(nbdof,numdof_global,true);

    if (verbose)
        std::cout<<"nbdof global "<<dist_index->getGlobalIndexSize()<<std::endl;

    return dist_index_shr_ptr;
}

//===== xDistStateDofCreator =======================================================================================================
template < typename Condition, typename VT >
xDistStateDofCreator <Condition, VT >::xDistStateDofCreator(xValueManagerDist<VT> & dm_, const std::string & sub) : dm(dm_), part_man(dm.getPartitionManager()),subset(sub) {}

//---------------------------------------
template < typename Condition, typename VT >
xDistStateDofCreator <Condition, VT >::xDistStateDofCreator(const typename xValueManagerDist<VT>::partman_t &part_man_,xValueManagerDist<VT> & dm_, const std::string & sub)
    : dm(dm_), part_man(part_man_),subset(sub) {}

//---------------------------------------
template < typename Condition, typename VT >
xValue < VT > * xDistStateDofCreator < Condition, VT >::operator()(const xValKey & key, xValue < VT >* v)
{
    // Usual condition : if state is already set we do nothing and return non null value v
    if (v->getState()) return v;
    // From Condition use test method to take a decision: return null if not to be treated
    if (!this->test(key)) return nullptr;
    // Extra condition related to distributed double manager: only key localy owned  are treated
    if (!part_man.getConstPartitionObject(key).isOwner()) return nullptr;

    // All condition are fulfill we can create dof state
    // Create new local numdof from subset size. Has only owned key are possible at this stage there is not problem
    // to use local numeration.
    int numdof = dm.size(subset) + 1;
    // Add value to subset
    dm.add(v, subset);
    // create state
    v->setState(new xStateOfValueDof(numdof));
    // Return non null value v
    return v;
}

//---------------------------------------
template < typename Condition, typename VT >
const typename xValueManagerDist<VT>::partman_t & xDistStateDofCreator < Condition, VT >::getPartitionManager()
{
    return part_man;
}

//---------------------------------------
template < typename Condition, typename VT >
xValueManagerDist<VT> &  xDistStateDofCreator < Condition, VT >::getValueManager()
{
    return dm;
}

//---------------------------------------
template < typename Condition, typename VT >
const std::string &  xDistStateDofCreator < Condition, VT >::getSubset()
{
    return subset;
}

//===== xDistStateDofCreatorUniform ====================================================================================================
template < typename PM, typename Condition, typename VT >
xDistStateDofCreatorUniform < PM, Condition, VT >::xDistStateDofCreatorUniform(const PM &part_man_, xValueManagerDist<VT> & val, const std::string & sub) : part_man(part_man_),val_manager(val), subset(sub) {}

//---------------------------------------
template < typename PM, typename Condition, typename VT >
xValue < VT > * xDistStateDofCreatorUniform < PM, Condition, VT >::operator()(const xValKey & key, xValue < VT >* v)
{
    // Usual condition : if state is already set we do nothing and return non null value v
    if (v->getState()) return v;
    // From Condition use test method to take a decision: return null if not to be treated
    if (!this->test(key)) return nullptr;
    // Extra condition related to distributed mesh: only key of local owned entity are treated
    AOMD::mEntity &e = *( key.getEnti());
    if (!part_man.getConstPartitionObject(e).isOwner()) return nullptr;

    // All condition are fulfill we can create dof state
    // Create new local numdof from subset size. Has only owned key are possible at this stage there is not problem
    // to use local numeration.
    int numdof = val_manager.size(subset) + 1;
    // Add value to subset
    val_manager.add(v, subset);
    // create state
    v->setState(new xStateOfValueDof(numdof));
    // Return non null value v
    return v;
}

//---------------------------------------
template < typename PM, typename Condition, typename VT >
const PM & xDistStateDofCreatorUniform < PM, Condition, VT >::getPartitionManager()
{
    return part_man;
}

//---------------------------------------
template < typename PM, typename Condition, typename VT >
xValueManagerDist<VT> &  xDistStateDofCreatorUniform < PM,Condition, VT >::getValueManager()
{
    return val_manager;
}

//---------------------------------------
template < typename PM, typename Condition, typename VT >
const std::string &  xDistStateDofCreatorUniform < PM,Condition, VT >::getSubset()
{
    return subset;
}

//===== xDistDofKeyAndInformationManagerUniform ==============================================================================================
template < typename PM, typename VT >
xDistDofKeyAndInformationManagerUniform < PM, VT >::xDistDofKeyAndInformationManagerUniform(const PM &part_man_, xValueManagerDist<VT> & dm_) : part_man(part_man_),dm(dm_) {}
//---------------------------------------
template < typename PM, typename VT >
xtool::xConstPartitionObject < AOMD::mEntity > xDistDofKeyAndInformationManagerUniform < PM, VT >::getConstPartitionObject(const data_t &o)
{
    AOMD::mEntity &e = *o.first.getEnti();
    return part_man.getConstPartitionObject(e);
}

//---------------------------------------
template < typename PM, typename VT >
auto xDistDofKeyAndInformationManagerUniform < PM, VT >::localObjectKey( const data_t &o)->information_key_t
{
    return o.first;
}

//---------------------------------------
template < typename PM, typename VT >
auto xDistDofKeyAndInformationManagerUniform < PM, VT >::remoteObjectKey(const xtool::xRemoteObject < AOMD::mEntity > &ro, const data_t &lo )->information_key_t
{
    xValKey remot_k = lo.first;
    remot_k.setEnti((AOMD::mEntity *) ( ro.getObjectAddress()));
    return remot_k;
}

//---------------------------------------
template < typename PM, typename VT >
auto xDistDofKeyAndInformationManagerUniform < PM, VT >::getInfo(const information_key_t &key, int sendto)->information_t
{
    xfem::xStateOfValue *state = ( dm.getValPtr(key) )->getState();
    assert(state);
    if (state->state == xfem::xStateOfValue::DOF)
    {
        xfem::xStateOfValueDof *dof = static_cast < xfem::xStateOfValueDof * >( state );
        return dof->Numdof;
    }
    else
        return -1;
}

//---------------------------------------
template < typename PM, typename VT >
void xDistDofKeyAndInformationManagerUniform < PM, VT >::setInfo(information_key_t key, const information_t &info, int receivedfrom)
{
    if (info > -1)
    {
        xfem::xStateOfValueDof*state = static_cast < xfem::xStateOfValueDof * >(( dm.getValPtr(key) )->getState());
        assert(state);
        state->Numdof = info;
    }
#ifndef NDEBUG
    else
    {
        xfem::xStateOfValue *state = ( dm.getValPtr(key) )->getState();
        assert(state);
        assert (state->state != xfem::xStateOfValue::DOF);
    }
#endif
}
//===== xDistFixedKeyManagerUniform =================================================================================================
template < typename PM, typename VT >
xDistFixedKeyManagerUniform < PM, VT >::xDistFixedKeyManagerUniform(const PM &part_man_, xValueManagerDist<VT> & dm_) : part_man(part_man_),dm(dm_) {}
//---------------------------------------
template < typename PM, typename VT >
auto xDistFixedKeyManagerUniform < PM, VT >::localObjectKey( const data_t &o)->information_key_t
{
    return o.first;
}

//---------------------------------------
template < typename PM, typename VT >
std::set < int >  xDistFixedKeyManagerUniform < PM, VT >::getMessageRanks( const information_key_t &key)
{
    // key is the xValkey to treat

    // first get entity of this key => local object
    AOMD::mEntity *lo = key.getEnti();

    // get remote object
    xtool::xConstPartitionObject < AOMD::mEntity > po = part_man.getConstPartitionObject(*lo);

    // int container (will be returned by this method), empty by default (no comunication)
    std::set < int > res;

    // if object have remote copy
    if (po.hasRemoteObject())
    {
        // we are now checking that it is a fixed value
        xfem::xStateOfValue *state = ( dm.getValPtr(key) )->getState();
        if (state && ( state->state == xfem::xStateOfValue::FIXED ))
        {
            // grabe remote proc id for object and store them in res
            for (const auto ro : po.getRemoteObjectsCollectionRange( ))
            {
                res.insert(ro.getProcId());
            }
        }
    }

    return res;

}
//===== xDistFixedInfoManagerUniform =================================================================================================
template < typename PM, typename VT >
xDistFixedInfoManagerUniform < PM, VT >::xDistFixedInfoManagerUniform(const PM &part_man_, xValueManagerDist<VT> & dm_) : part_man(part_man_),dm(dm_),nbc(0) {}
//---------------------------------------
template < typename PM, typename VT >
void xDistFixedInfoManagerUniform < PM, VT >::getInfo(information_key_t key, xtool::xMpiInputBuffer & buff, int sendto)
{
    // associated local object
    AOMD::mEntity *lo = key.getEnti();

    // get remote object
    xtool::xConstPartitionObject < AOMD::mEntity > po = part_man.getConstPartitionObject(*lo);
    assert (po.hasRemoteObject());
    const AOMD::mEntity *ro = po.getRemoteObjectOn(sendto);
    //pack the remote object
    buff.pack(&ro, 1, MPI_AINT);

    // pack now the rest of the key !!!
    // here for this first implementation we choose to transfers string. Less expensive would have been to use
    // associated xValKey::ids_size_t (phys/geom) but we don't have for the moment the insurance that those information are always created
    // in the same order in all proc. If someone visit xtool::xStringManager so that it work in // nicely then use of  xValKey::ids_size_t will be fine.
    // Experimental version with use of hash function is a promising solution which may lead to use xValKey::ids_size_t
    // the geom string
    std::string geom = xKeyInfo::getGeomName(key.getGeom());
    buff.pack(geom.c_str(), geom.size()+1, MPI_CHAR);
    // the phys string
    std::string phys = xKeyInfo::getPhysName(key.getPhys());
    buff.pack(phys.c_str(), phys.size()+1, MPI_CHAR);

    // Last but not the least, the value fixed
    VT v = ( dm.getValPtr(key) )->getVal();
    buff.pack(&v, 1, MPI_DOUBLE);
}
//---------------------------------------
template < typename PM, typename VT >
void xDistFixedInfoManagerUniform < PM, VT >::setInfo(const xtool::xMpiOutputBuffer & buff, int receivedfrom)
{
    char c;
    VT v;
    AOMD::mEntity *lo;
    while (!buff.exhausted() )
    {
        // get entity
        buff.unPack(&lo, 1, MPI_AINT);

        // get geom and phys
        std::string geom,phys;
        buff.unPack(&c, 1, MPI_CHAR);
        while (c != '\0')
        {
            geom += c;
            buff.unPack(&c, 1, MPI_CHAR);
        }
        buff.unPack(&c, 1, MPI_CHAR);
        while (c != '\0')
        {
            phys += c;
            buff.unPack(&c, 1, MPI_CHAR);
        }

        // create key
        xValKey key(xKeyInfo::getPhysId(phys),xKeyInfo::getGeomId(geom),lo);

        // get val
        auto val = dm.getValPtr(key);
        buff.unPack(&v, 1, MPI_DOUBLE);

        // check state of value
        xfem::xStateOfValue *state = val->getState();
        if (!state)
        {
            val->setVal(v);
            val->setState(new xStateOfValueFixed<VT>(val));
            ++nbc;
        }
        else if (state->state != xfem::xStateOfValue::FIXED)
        {
            throw -1;
        }
    }
}
//---------------------------------------
template < typename PM, typename VT >
size_t xDistFixedInfoManagerUniform < PM, VT >::getApproxDataSize() const
{
    // taking this example :
    // DISPLACEMENT_Y Lagrange_Order_1_Node
    return ( sizeof( VT ) + sizeof( AOMD::mEntity * ) + sizeof( char )*37 );
}




//===== xDistDofKeyAndInformationManager =================================================================================================
template<typename VT>
xfem::xDistDofKeyAndInformationManager<VT>::xDistDofKeyAndInformationManager(xValueManagerDist<VT> & dm_,const typename xValueManagerDist<VT>::partman_t &pm)
    : dm(dm_),part_man(pm) {}
//---------------------------------------

template<typename VT>
xtool::xConstPartitionObject < xfem::xValKey > xfem::xDistDofKeyAndInformationManager<VT>::getConstPartitionObject(const data_t &o)
{
    return part_man.getConstPartitionObject(o.first);
}

//---------------------------------------

template<typename VT>
auto xfem::xDistDofKeyAndInformationManager<VT>::localObjectKey( const data_t &o)->information_key_t
{
    const xValKey *lk = dm.findKeyAdress(o.first);
    assert(lk);
    return lk;
}

//---------------------------------------
template<typename VT>
auto xfem::xDistDofKeyAndInformationManager<VT>::remoteObjectKey(const xtool::xRemoteObject < xValKey > &ro, const data_t &lo )->information_key_t
{
    return static_cast < information_key_t >( ro.getObjectAddress());
}

//---------------------------------------

template<typename VT>
auto xfem::xDistDofKeyAndInformationManager<VT>::getInfo(const information_key_t &key, int sendto)->information_t
{
    xfem::xStateOfValue *state = ( dm.getValPtr(*key) )->getState();
    assert(state);
    if (state->state == xfem::xStateOfValue::DOF)
    {
        xfem::xStateOfValueDof *dof = static_cast < xfem::xStateOfValueDof * >( state );
        return dof->Numdof;
    }
    else
        return -1;
}

//---------------------------------------

template<typename VT>
void xfem::xDistDofKeyAndInformationManager<VT>::setInfo(information_key_t key, const information_t &info, int receivedfrom)
{
    if (info > -1)
    {
        xfem::xStateOfValueDof*state = static_cast < xfem::xStateOfValueDof * >(( dm.getValPtr(*key) )->getState());
        assert(state);
        state->Numdof = info;
    }
#ifndef NDEBUG
    else
    {
        xfem::xStateOfValue *state = ( dm.getValPtr(*key) )->getState();
        assert(state);
        assert (state->state != xfem::xStateOfValue::DOF);
    }
#endif
}
//===== xDistFixedKeyManager =================================================================================================
template <typename VT>
xfem::xDistFixedKeyManager<VT>::xDistFixedKeyManager(xValueManagerDist<VT> & dm_) : dm(dm_),part_man(dm.getPartitionManager()) {}
//---------------------------------------
template <typename VT>
auto xfem::xDistFixedKeyManager<VT>::localObjectKey( const data_t &o)->information_key_t
{
    return(o.first);
}

//---------------------------------------
template <typename VT>
std::set < int >  xfem::xDistFixedKeyManager<VT>::getMessageRanks( const information_key_t &key)
{
    // key is the xValkey to treat

    // get remote object
    xtool::xConstPartitionObject < xfem::xValKey > po = part_man.getConstPartitionObject(key);

    // int container (will be returned by this method), empty by default (no comunication)
    std::set < int > res;

    // if object have remote copy
    if (po.hasRemoteObject())
    {
        // we are now checking that it is a fixed value
        xfem::xStateOfValue *state = ( dm.getValPtr(key) )->getState();
        if (state && ( state->state == xfem::xStateOfValue::FIXED ))
        {
            // grabe remote proc id for object and store them in res
            for (const auto ro : po.getRemoteObjectsCollectionRange( ))
            {
                res.insert(ro.getProcId());
            }
        }
    }

    return res;

}
//===== xDistFixedInfoManager =================================================================================================
template <typename VT>
xfem::xDistFixedInfoManager<VT>::xDistFixedInfoManager(xValueManagerDist<VT> & dm_) : dm(dm_),part_man(dm.getPartitionManager()),nbc(0) {}
//---------------------------------------
template <typename VT>
auto  xfem::xDistFixedInfoManager<VT>::getInfo(information_key_t key, int sendto) -> information_t
{
    // Get remote object
    xtool::xConstPartitionObject < xValKey > po = part_man.getConstPartitionObject(key);
    assert (po.hasRemoteObject());

    // information to send
    information_t info;

    // The remote xValKey adress
    info.first= po.getRemoteObjectOn(sendto);

    // The fixed value
    info.second= ( dm.getValPtr(key) )->getVal();

    return info;
}
//---------------------------------------
template <typename VT>
void xfem::xDistFixedInfoManager<VT>::setInfo( const std::vector < information_t > &infos, int receivedfrom)
{
    for (auto info : infos)
    {
        // get val
        auto val = dm.getValPtr(*info.first);

        // check state of value
        xfem::xStateOfValue *state = val->getState();
        if (!state)
        {
            val->setVal(info.second);
            val->setState(new xStateOfValueFixed<VT>(val));
            ++nbc;
        }
        else if (state->state != xfem::xStateOfValue::FIXED)
        {
            throw -1;
        }
    }
}






}         //end namespace
#endif
