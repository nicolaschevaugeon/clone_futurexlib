/*
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
 */
#ifdef _XDISTBLOCKDENSEMATRIX_H
#ifndef _XDISTBLOCKDENSEMATRIX_IMP_H
#define _XDISTBLOCKDENSEMATRIX_IMP_H

namespace xlinalg
{
//===== xKeyManagerDistBlockDenseMatrix =====================================================================================
template < typename T, typename DEFINED, typename PATTERN >
xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xKeyManagerDistBlockDenseMatrix::xKeyManagerDistBlockDenseMatrix (const xDistIndex &dist_index_,xDistBlockDenseMatrix::Vector &data_) :
    dist_index(dist_index_)
    ,dat(data_)
{}
template < typename T, typename DEFINED, typename PATTERN >
auto xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xKeyManagerDistBlockDenseMatrix::localObjectKey( const data_t &lo)->information_key_t
{
    return &dat[dist_index.getPackedIndex(lo)];
}
template < typename T, typename DEFINED, typename PATTERN >
xtool::xConstPartitionObject < xDistIndex::idx_t > xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xKeyManagerDistBlockDenseMatrix::getConstPartitionObject( const data_t &lo)
{
    return dist_index.getConstPartitionObject(lo);
}
template < typename T, typename DEFINED, typename PATTERN >
auto xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xKeyManagerDistBlockDenseMatrix::remoteObjectKey(const xtool::xRemoteObject < xDistIndex::idx_t > &ro, const data_t &lo)->information_key_t
{

    // Here we cheat : normally we should have given remote address of the first column cell corresponding to lo like in localObjectKey
    // we give a proper pointer address. But this would have required to create a partition manager just for this instance vector
    // and exchange address. Not pleasant at all.
    // In fact all what we need here is something giving the right order in remote proc of impacted cell in dat .....
    // For that we use address of remote index of dist_index instance corresponding to lo given by ro. Has those address represent
    // pack_id address these are well ordered. Remote will use for keys it's local vector address constructed by use of getPackedIndex,
    // says the pack_id address ....
    //
    // To make this remote address a key we need to cast it.
    return reinterpret_cast < information_key_t >( const_cast < xDistIndex::idx_t * >( ro.getObjectAddress()));
}
//===== xInfoManagerDistBlockDenseMatrixReduce  =====================================================================================
template < typename T, typename DEFINED, typename PATTERN >
xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xInfoManagerDistBlockDenseMatrixReduce::xInfoManagerDistBlockDenseMatrixReduce(Vector &dat_,index_column_asso_container_t &exchanged_column_terms_) :
    dat(dat_)
    ,exchanged_column_terms(exchanged_column_terms_)
    ,zero(xtool::xDataType < T >::zero())
    ,mpi_data_t(xtool::xMPIDataType < T >( ))
    ,moy(0)
{}
template < typename T, typename DEFINED, typename PATTERN >
void xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xInfoManagerDistBlockDenseMatrixReduce::getInfo(information_key_t key, xtool::xMpiInputBuffer & buff,int sendto)
{
    // in reduce key will allways corespond to "not owned" lines that we want to add in processor owning it
    const size_t i = key-&dat[0];
    auto itf = exchanged_column_terms.find(std::make_pair(i,sendto));
    if (itf != exchanged_column_terms.end())
    {
        for (auto p : itf->second)
        {
            T &val = *( p.second );
            buff.pack(p.second,1,mpi_data_t);
            val = zero;
        }
    }
    else
    {
        throw -1;
    }
}
template < typename T, typename DEFINED, typename PATTERN >
void xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xInfoManagerDistBlockDenseMatrixReduce::setInfo(information_key_t key,  const xtool::xMpiOutputBuffer & buff, int receivedfrom)
{
    const size_t i = key-&dat[0];
    auto itf = exchanged_column_terms.find(std::make_pair(i,receivedfrom));
    if (itf != exchanged_column_terms.end())
    {
        for (auto p : itf->second)
        {
            T val;
            buff.unPack(&val,1,mpi_data_t);
            *( p.second ) += val;
        }
    }
    else
    {
        throw -1;
    }
}
template < typename T, typename DEFINED, typename PATTERN >
void xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xInfoManagerDistBlockDenseMatrixReduce::setMoy(int moy_)
{
    moy = moy_;
}
template < typename T, typename DEFINED, typename PATTERN >
size_t xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xInfoManagerDistBlockDenseMatrixReduce::getApproxDataSize(void)
{
    return sizeof( double )*moy;
}
//===== xInfoManagerDistBlockDenseMatrixSet  ========================================================================================
template < typename T, typename DEFINED, typename PATTERN >
xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xInfoManagerDistBlockDenseMatrixSet::xInfoManagerDistBlockDenseMatrixSet(Vector &dat_,index_column_asso_container_t &exchanged_column_terms_) :
    dat(dat_)
    ,exchanged_column_terms(exchanged_column_terms_)
    ,mpi_data_t(xtool::xMPIDataType < T >( ))
    ,moy(0)
{}
template < typename T, typename DEFINED, typename PATTERN >
void xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xInfoManagerDistBlockDenseMatrixSet::getInfo(information_key_t key, xtool::xMpiInputBuffer & buff,int sendto)
{
    const size_t i = key-&dat[0];
    auto itf = exchanged_column_terms.find(std::make_pair(i,sendto));
    if (itf != exchanged_column_terms.end())
    {
        for (auto p : itf->second) buff.pack(p.second,1,mpi_data_t);
    }
    else
    {
        throw -1;
    }
}
template < typename T, typename DEFINED, typename PATTERN >
void xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xInfoManagerDistBlockDenseMatrixSet::setInfo(information_key_t key,  const xtool::xMpiOutputBuffer & buff, int receivedfrom)
{
    const size_t i = key-&dat[0];
    auto itf = exchanged_column_terms.find(std::make_pair(i,receivedfrom));
    if (itf != exchanged_column_terms.end())
    {
        for (auto p : itf->second) buff.unPack(p.second,1,mpi_data_t);
    }
    else
    {
        throw -1;
    }
}
template < typename T, typename DEFINED, typename PATTERN >
void xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xInfoManagerDistBlockDenseMatrixSet::setMoy(int moy_)
{
    moy = moy_;
}
template < typename T, typename DEFINED, typename PATTERN >
size_t xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xInfoManagerDistBlockDenseMatrixSet::getApproxDataSize(void)
{
    return sizeof( double )*moy;
}
//===== xDistBlockDenseMatrix ===========================================================================================
template < typename T, typename DEFINED, typename PATTERN >
xDistBlockDenseMatrix < T,DEFINED,PATTERN >::xDistBlockDenseMatrix( const xDistIndex & dist_index_ ) :
    n_bloc(dist_index_.getPackedIndexSize())
    ,dat(n_bloc*n_bloc, xtool::xDataType < T >::zero())
    ,dist_index(dist_index_)
    ,reduce_info(dat,exchanged_column_terms)
    ,set_info(dat,exchanged_column_terms)
    ,reduce_keys(dist_index.getComm() )
    ,set_keys(dist_index.getComm() )
    ,status(RESET)
    ,proc_id(-1)
{

    // Reduce operation consite in adding "not owned" lines into processor that own it.
    // From this line only "not owned" column have to be send for the proc which own the line.
    //
    // Set operation consite in setting "owned" lines into processor that do no own it but have part of it.
    // From this line only "owned" column have to be send for the proc which do not own the line.
    //
    // In consequence tables must be set to give column index per line per remote proc for this two situation.
    // No comunication is needed as we use dist_index which is giving correct remote order

    // loop on line to fill exchanged_column_terms need by both reduce_info and set_info
    for (auto idi : dist_index)
    {
        const auto &poi = dist_index.getConstPartitionObject(idi);
        const int i = dist_index.getPackedIndex(idi);
        // We are only interested by line with remote in some proc
        if (poi.hasRemoteObject())
        {
            // For Set: we are getting only lines that are owned
            // For Reduce: only lines that are owned are modified
            if (poi.isOwner())
            {
                for (auto roi : poi.getRemoteObjectsCollectionRange())
                {
                    int remot_proc_i = roi.getProcId();
                    auto key = std::make_pair(i,remot_proc_i);

                    for (auto idj : dist_index)
                    {
                        const auto &poj = dist_index.getConstPartitionObject(idj);
                        // only owned column interest us
                        if (poj.isOwner())
                        {
                            for (auto roj : poj.getRemoteObjectsCollectionRange())
                            {
                                // only column with a counterpart to line remote proc interest us
                                if (roj.getProcId() == remot_proc_i)
                                {
                                    T *p = &dat[i+n_bloc*dist_index.getPackedIndex(idj)];
                                    auto itf = exchanged_column_terms.find(key);
                                    if ( itf != exchanged_column_terms.end())
                                    {
                                        itf->second.insert(std::make_pair(dist_index.getPackedAdress(idj),p));
                                    }
                                    else
                                    {
                                        index_column_asso_t tmp;
                                        tmp.insert ( std::make_pair(dist_index.getPackedAdress(idj),p));
                                        exchanged_column_terms.insert(std::make_pair(key,tmp));
                                    }
                                    break;
                                }
                            }

                        }
                    }
                }
            }
            // For Set: we are setting only lines that are not owned
            // For Reduce: we are sending only lines that are not owned
            else
            {
                int remot_proc_i_owner = poi.getOwner().getProcId();
                auto key = std::make_pair(i,remot_proc_i_owner);
                for (auto idj : dist_index)
                {
                    const auto &poj = dist_index.getConstPartitionObject(idj);
                    // only not owned column interest us
                    if (!poj.isOwner())
                    {
                        auto rojo = poj.getOwner();
                        int remote_proc_j_owner = rojo.getProcId();
                        // only column which correspond to proc owning line and which is also owner itself of the column interest us, if any
                        if (remote_proc_j_owner == remot_proc_i_owner )
                        {
                            T *p = &dat[i+n_bloc*dist_index.getPackedIndex(idj)];
                            auto itf = exchanged_column_terms.find(key);
                            if ( itf != exchanged_column_terms.end())
                            {
                                itf->second.insert(std::make_pair(static_cast < const idx_t * >( rojo.getObjectAddress()),p));
                            }
                            else
                            {
                                index_column_asso_t tmp;
                                tmp.insert(std::make_pair(static_cast < const idx_t * >( rojo.getObjectAddress()),p));
                                exchanged_column_terms.insert(std::make_pair(key,tmp));
                            }
                        }
                    }

                }
            }
        }
    }
    int ms = 0;
    for (auto asso : exchanged_column_terms)
    {
        ms += asso.second.size();
    }
    int moy = ms/exchanged_column_terms.size();
    reduce_info.setMoy(moy);
    set_info.setMoy(moy);
//create key manager
    xKeyManagerDistBlockDenseMatrix key_manager(dist_index,dat);
//create keys
    reduce_keys.accumulateKeysOwnerGather(dist_index.begin(),dist_index.end(),key_manager);
    set_keys.accumulateKeysOwnerScatter(dist_index.begin(),dist_index.end(),key_manager);
// getting rank in world
    MPI_Comm_rank(dist_index.getComm(),&proc_id);
}
template < typename T, typename DEFINED, typename PATTERN >
inline void xDistBlockDenseMatrix < T,DEFINED,PATTERN >::Printdata( std::ostream  & out) const
{
    out << "###########################"<<std::endl;
    out << "Matrix of dimension " << dist_index.getGlobalIndexSize() << std::endl;
    out << "Local contribution to this matrix is:\n";
    for (auto idi : dist_index )
    {
        out<<" | "<<std::setprecision(5)<<std::setw(5)<<idi;
    }
    out<<" | "<<std::endl;
    for (int i = 0; i < n_bloc; ++i) out<<std::setw(8)<<"-";
    out<<std::endl;
    for (int i = 0; i < n_bloc; ++i)
    {
        for (int j = 0; j < n_bloc; ++j)
        {
            out<<" | "<<std::setw(5)<<dat[i+n_bloc*j];
        }
        out<<" | "<<std::endl;
    }
    out << std::endl;
    out << "###########################"<<std::endl;
}
template < typename T, typename DEFINED, typename PATTERN >
void xDistBlockDenseMatrix < T,DEFINED,PATTERN >::reduceOnOwner()
{
    if (!( status&REDUCED ))
    {
        xtool::exchangeInformation(reduce_keys,reduce_info);
    }
    status |= REDUCED;
}

template < typename T, typename DEFINED, typename PATTERN >
void xDistBlockDenseMatrix < T,DEFINED,PATTERN >::getMemoryAccess(int &ng,int &nl,T **data)
{
    ng = dist_index.getGlobalIndexSize();
    nl = n_bloc;
    *data = &dat[0];
}
template < typename T, typename DEFINED, typename PATTERN >
void xDistBlockDenseMatrix < T,DEFINED,PATTERN >::gemv( T alpha, T beta,const xlinalg::xDistVector<T> & X, xlinalg::xDistVector<T> & Y)
{
    // switch X in insert mode off
    bool X_imon = X.isInInsertModeOn();
    if (X_imon)
        const_cast < xlinalg::xDistVector<T> & >( X ).switchInsertModeOff();

    // swith to global state so that we can use term correctely
    const_cast < xlinalg::xDistVector<T> & >( X ).switchToGlobalValue();

    // switch Y in insert mode on
    bool Y_imon = Y.isInInsertModeOn();
    if (!Y_imon)
        Y.switchInsertModeOn();

    // compute the matrix vector product
    //if ( status&REDUCED )
    //{
    // TODO: if matrix is reduced some none owned diagonal block are null, skip their computation.
    //}
    //else
    xPolicyxDistBlockDenseMatrixGemv < T,PATTERN >::mv(n_bloc, alpha, &dat[0], n_bloc, &X[0], beta, &Y[0] );

    // revert to initial state
    if (!Y_imon)
        Y.switchInsertModeOff();
    if (X_imon)
        const_cast < xlinalg::xDistVector<T> & >( X ).switchInsertModeOn();
}


//===== xPolicyxDistBlockDenseMatrixGemv ================================================================================
template < >
class xPolicyxDistBlockDenseMatrixGemv < double,xTraitMatrixLowerSym >
{
    public:
        static void   mv( const int & n, const double & alpha, const double*A, const int & LDA,const double *X, const double & beta, double *Y);
};
template < >
class xPolicyxDistBlockDenseMatrixGemv < double,xTraitMatrixUpperSym >
{
    public:
        static void   mv( const int & n, const double & alpha, const double*A, const int & LDA,const double *X, const double & beta, double *Y);
};
template < >
class xPolicyxDistBlockDenseMatrixGemv < double,xTraitMatrixUnSym >
{
    public:
        static void   mv( const int & n, const double & alpha, const double*A, const int & LDA,const double *X, const double & beta, double *Y);
};

}                                                      // end of namespace
#endif
#endif







