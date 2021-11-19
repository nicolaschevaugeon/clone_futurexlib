/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifdef _xPartitionManagerTools_
#ifndef _xPartitionManagerTools_imp_
#define _xPartitionManagerTools_imp_
namespace xtool
{
/// class specifique to createPartitionManagerForSubGroup function
namespace createPartitionManagerForSubGroupFunc
{
template < typename PM, typename T >
class InfoManagerFilterSGPMBase
{
    public:
        /// basic constructor which retrive partition manager associated to this class
        InfoManagerFilterSGPMBase( const PM &partman_ ) : partman(partman_){}

        // mandatory types for information class family
        typedef xtool::homogeneous_data_style_trait data_style_trait;
        typedef xtool::send_only_keys_communication_trait communication_trait;
        typedef const T * information_key_t;

        // mandary method (partial) for information class familly
        information_key_t localObjectKey( const T &o)
        {
            return &o;
        }
        information_key_t remoteObjectKey(const xtool::xRemoteObject < T > &ro, const T &lo )
        {
            return ro.getObjectAddress();
        }
        xtool::xConstPartitionObject < T > getConstPartitionObject(const T &e)
        {
            return partman.getConstPartitionObject(e);
        }

    protected:
        const PM &partman;
};
//  filter partition manager : needed after first local work
//! Filtration is done in 2 pass => 2 info manager => one feeding the other one
template < typename PM, typename T, template < typename > class DM >
class InfoManagerFilterSGPM1 : public InfoManagerFilterSGPMBase < PM,T >
{
    public:
        typedef DM < std::set < std::pair < int,const T * > > > data_man_t;
        typedef typename InfoManagerFilterSGPMBase < PM,T >::information_key_t information_key_t;
        typedef std::pair < const T *,const T * > information_t;

        /// Constructor
        InfoManagerFilterSGPM1(const PM &sg_part_man, data_man_t &to_destroy_) :
            InfoManagerFilterSGPMBase < PM,T >(sg_part_man)
            ,to_destroy(to_destroy_)
        {}



        // last mandatory method for information class familly
        information_t getInfo(information_key_t key, int sendto)
        {
            // lk is the the key : pointer to local entity object
            // information is the remote object coresponding to sendto
            // get partion object
            auto po = ( const_cast < PM & >( InfoManagerFilterSGPMBase < PM,T >::partman )).getConstPartitionObject(*key);
            assert (po.hasRemoteObject());

            // get remote object
            auto ro = po.getRemoteObjectOn(sendto);
            assert(ro);

            //return std::make_pair(const_cast < information_t >( ro ),lk);
            return std::make_pair(ro,key);

        }
        void setInfo(const std::vector < information_t > &infos, int receivedfrom)
        {
            for (auto i : infos)
            {
                // xConstPartitionObject object related to this entity present in remote proc receivedfrom
                auto po = ( const_cast < PM & >( InfoManagerFilterSGPMBase < PM,T >::partman )).getConstPartitionObject(*( i.first ));
                // This entity do not have any counterpart at sg level but receive from another proc. It must be
                // inserted in to destroy so that reverse communication remove the received remote connection
                if (!po.getRemoteObjectOn(receivedfrom))
                {
                    to_destroy.setData(*( const_cast < T * >( i.first ))).insert(std::make_pair(receivedfrom,i.second));
                }
            }

        }
        std::set < int > getMessageRanks( information_key_t lk)
        {
            // lk is the pointer to local object
            // So we can used it to access to remote object
            auto po = ( const_cast < PM & >( InfoManagerFilterSGPMBase < PM,T >::partman )).getConstPartitionObject(*lk);
            assert (po.hasRemoteObject());

            // grabe remote proc id for object and store them in res
            std::set < int > res;
            for (const auto ro : po.getRemoteObjectsCollectionRange( )) res.insert(ro.getProcId());

            return res;
        }
    private:
        data_man_t &to_destroy;
};
template < typename PM, typename T, template < typename > class DM >
class InfoManagerFilterSGPM2 : public InfoManagerFilterSGPMBase < PM,T >
{
    public:
        /// Constructor
        InfoManagerFilterSGPM2(const PM &sg_part_man, typename InfoManagerFilterSGPM1 < PM,T,DM >::data_man_t &to_destroy_) :
            InfoManagerFilterSGPMBase < PM,T >(sg_part_man)
            ,to_destroy(to_destroy_)
        {}

        // last mandatory types for information class familly
        typedef const T * information_t;
        typedef typename InfoManagerFilterSGPMBase < PM,T >::information_key_t information_key_t;

        // last mandatory method for information class familly
        information_t getInfo(information_key_t key, int sendto)
        {
            assert(to_destroy.getData(*key));
            for ( auto i : *( to_destroy.getData(*key)) )
                if (i.first == sendto) return i.second;
            throw -1;
        }
        void setInfo(const std::vector < information_t > &infos, int receivedfrom)
        {
            for (auto i : infos)
            {
                // xPartitionObject object related to this entity which must be cleaned
                auto po = ( const_cast < PM & >( InfoManagerFilterSGPMBase < PM,T >::partman )).getPartitionObject(*const_cast < T * >( i ));
                assert (po.getRemoteObjectOn(receivedfrom));
                po.remove(receivedfrom);
            }

        }
        std::set < int > getMessageRanks( information_key_t lk)
        {
            // lk is the pointer to local object
            // check on non existence of remote copy
            assert (!( const_cast < PM & >( InfoManagerFilterSGPMBase < PM,T >::partman )).getConstPartitionObject(*lk).hasRemoteObject());

            std::set < int > res;
            // get set of int where to send a destroy information : a remote object address to be cleaned for this proc id
            assert(to_destroy.getData(*lk));
            for ( auto i : *( to_destroy.getData(*lk)) ) res.insert(i.first);

            return res;
        }
    private:
        typename InfoManagerFilterSGPM1 < PM,T,DM >::data_man_t &to_destroy;
};
};
// Implementation
template <  template < typename > class DM, typename T, typename RANGE, typename GPM, typename SGPM >
void createPartitionManagerForSubGroup(RANGE sub_group_boundary_range, const GPM & general_partman, SGPM & sub_group_partman)
{
    // communicator : general and sub group
    MPI_Comm univ = general_partman.getComm();
    MPI_Comm sg_univ = sub_group_partman.getComm();
    assert(univ != MPI_COMM_NULL);
    assert(sg_univ != MPI_COMM_NULL);

    int nb_proc;
    MPI_Comm_size(univ,&nb_proc);

    // If underlying mesh is distributed (nb_proc>1), sub group may be on more then one proc
    // Otherwise it is finished as the group will not comunicate with any one else and his partition manager has
    // to remain empty
    // Note: with some comunicator sub group may exist in different process but may be treated as local group for those
    // process because each process do have a comuniquator with a unique prosses  (MPI_COMM_SELF is an example
    // of such situation). This case must be treated the same way as the other has comunication pattern involve all
    // proc of this comunicator.
    if (nb_proc > 1)
    {

        int same;
        MPI_Comm_compare(univ,sg_univ,&same);
        DM < bool > done;

        // if we deal with different comunicator (i.e. not having same rank) we have to filter with ranks
        if ( same != MPI_IDENT || same != MPI_CONGRUENT)
        {

            // group to set translation
            MPI_Group ggu, gu;
            MPI_Comm_group(univ, &gu );
            MPI_Comm_group(sg_univ, &ggu );


            // translated and none translated id
            std::vector < int > ranks(nb_proc);
            for (int i = 0; i < nb_proc; i++) ranks[i] = i;
            std::vector < int > sg_ranks(nb_proc,-1);
            MPI_Group_translate_ranks( gu, nb_proc, &ranks[0], ggu, &sg_ranks[0] );


            // loop on proc boundary of the group (given by user as a range) and create all remote information only if those
            // remote information are in the sg_univ and local info in group. This is isolating future communication only on entity of this group
            // To avoid unnecessary repetition use of a data manger to mark done entity

            for (auto eb : sub_group_boundary_range)
            {

                const int dimm1 = eb->getLevel();
                for (int i = 0; i < dimm1; ++i)
                {
                    for (int j = 0,m = eb->size(i); j < m; ++j)
                    {
                        T & ei = *( eb->get(i,j));
                        if (!done.getData(ei))
                        {
                            done.setData(ei) = true;
                            auto po = general_partman.getConstPartitionObject(ei);
                            auto sg_po = sub_group_partman.getPartitionObject(ei);
                            for (auto ro : po.getRemoteObjectsCollectionRange())
                            {
                                const int sg_proc_id = sg_ranks[ro.getProcId()];
                                if (sg_proc_id != MPI_UNDEFINED)
                                    sg_po.insert(sg_proc_id, ro.getObjectAddress());
                            }
                        }
                    }
                }
                auto po = general_partman.getConstPartitionObject(*eb);
                auto sg_po = sub_group_partman.getPartitionObject(*eb);
                for (auto ro : po.getRemoteObjectsCollectionRange())
                {
                    const int sg_proc_id = sg_ranks[ro.getProcId()];
                    if (sg_proc_id != MPI_UNDEFINED)
                        sg_po.insert(sg_proc_id, ro.getObjectAddress());
                }
            }

        }
        else
        {
            // loop on proc boundary of the group (given by user as a range) and create all remote information only if those
            // local information are in the group. This is isolating future communication only on entity of this group
            // To avoid unnecessary repetition use of a data manger to mark done entity
            for (auto eb : sub_group_boundary_range)
            {

                const int dimm1 = eb->getLevel();
                for (int i = 0; i < dimm1; ++i)
                {
                    for (int j = 0,m = eb->size(i); j < m; ++j)
                    {
                        auto & ei = *( eb->get(i,j));
                        if (!done.getData(ei))
                        {
                            done.setData(ei) = true;
                            auto po = general_partman.getConstPartitionObject(ei);
                            auto sg_po = sub_group_partman.getPartitionObject(ei);
                            for (auto ro : po.getRemoteObjectsCollectionRange())
                            {
                                sg_po.insert(ro.getProcId(), ro.getObjectAddress());
                            }
                        }
                    }
                }
                auto po = general_partman.getConstPartitionObject(*eb);
                auto sg_po = sub_group_partman.getPartitionObject(*eb);
                for (auto ro : po.getRemoteObjectsCollectionRange())
                {
                    sg_po.insert(ro.getProcId(), ro.getObjectAddress());
                }
            }
        }
        // Note : The selecting above algorithm is in default at least in the following case :
        //  For a given group
        //     a node A of sub_group_boundary_range in proc k is duplicate on i and j proc in sg_univ communicator
        //     this node is in sub_group_boundary_range of i proc
        //     this node is NOT in sub_group_boundary_range of j proc. Says it is in partition manager related to k proc but not in sub group ....
        //  in this case A in k will have 2 remotes (i,j)
        //               A in i will have 2 remotes (k,j)
        //               A in j will have no remote
        // which is incorrect and explain why an extra exchange is mandatory to remove extra relation in  partition manager for this  sub group
        // In example above refitting lead to
        //               A in k will have 1 remote (i)
        //               A in i will have 1 remote (k)
        //               A in j will have no remote
        // The first pass have done a huge reduction from general mesh on all proc to sub group on a possibly restricted mpi communicator.
        // Here we exchange only on sub group boundary entity. It is, a priori, a smaller set then general mesh and only
        // in between proc of the group which is may be a smaller set then general communicator.
        typename createPartitionManagerForSubGroupFunc::InfoManagerFilterSGPM1 < SGPM,T,DM >::data_man_t to_destroy;
        createPartitionManagerForSubGroupFunc::InfoManagerFilterSGPM1 < SGPM,T,DM > update1_sg_partman_key_info(sub_group_partman,to_destroy);
        createPartitionManagerForSubGroupFunc::InfoManagerFilterSGPM2 < SGPM,T,DM > update2_sg_partman_key_info(sub_group_partman,to_destroy);
        xtool::xKeyContainerSendOrRecv < typename createPartitionManagerForSubGroupFunc::InfoManagerFilterSGPM1 < SGPM,T,DM >::information_key_t >  update_sg_partman_key_container(sg_univ);

        // first pass : all boundary proc sub group are exchanged in a key send only fashion => identification of remote
        // present only in one side
        update_sg_partman_key_container.accumulateKeys(sub_group_partman.beginObject(),sub_group_partman.endObject(),update1_sg_partman_key_info);
        xtool::exchangeInformation(update_sg_partman_key_container,update1_sg_partman_key_info);

        // Second pass : only key identified to be destroyed from sender (above pass) must be set and send back to
        // be destroy in sender
        // This set may be empty and is supposed to be small
        update_sg_partman_key_container.clearKeys();
        update_sg_partman_key_container.accumulateKeys(to_destroy.beginKey(),to_destroy.endKey(),update2_sg_partman_key_info);
        xtool::exchangeInformation(update_sg_partman_key_container,update2_sg_partman_key_info);
        /*
           // DEBUG DEBUG
           std::cout << "List of partman for patch id "<<id<<std::endl;
           for (auto ite = part_man_patch->beginObject(), itee = part_man_patch->endObject(); ite != itee; ++ite)
           {
            printEntity(*ite,false);
            std::cout << "Remote :";
            auto po = part_man_patch->getConstPartitionObject(**ite);
            if (po.hasRemoteObject())
            {
                for (auto ro : po.getRemoteObjectsCollectionRange( ))
                {
                    std::cout << " ("<<ro.getProcId()<<","<<ro.getObjectAddress()<<")";
                }
            }
            else
            {
                std::cout << " strange no remote";
            }

            std::cout << std::endl;

           }
         */

    }

}
template <  typename PM >
void exportInGraphvizDof(const PM & partman, std::string name, size_t filter)
{

    int proc_id,nb_proc;
    MPI_Comm world = partman.getComm();
    MPI_Comm_rank(world,&proc_id);
    MPI_Comm_size(world,&nb_proc);
    auto r = make_range(partman.beginObject(),partman.endObject());
    std::vector < int > sizefullmsg(proc_id ? 0 : nb_proc,0);
    std::vector < int > dispfullmsg(proc_id ? 0 : nb_proc,0);
    std::vector < char > msg;
    for (int dim = 0; dim < 3; ++dim)
    {
        msg.clear();
        for (auto e : r)
        {
            if (e->getLevel() == dim)
            {
                auto po = partman.getConstPartitionObject(*e);
                const auto & range = po.getRemoteObjectsCollectionRange();
                if (range.size() > filter)
                {
                    std::stringstream ss;
                    ss<<"p"<<proc_id<<"_"<<e<<" -> {";
                    int i = -1;
                    for (const auto ro : range)
                    {
                        if (++i) ss<<" ";
                        ss<<"p"<<ro.getProcId()<<"_"<<ro.getObjectAddress();
                    }
                    ss<<"}\n";
                    std::string msgi = ss.str();
                    msg.insert(msg.end(),msgi.begin(),msgi.end());
                }
            }
        }
        int size_msg = msg.size();
        MPI_Gather( &size_msg,1,MPI_INT,sizefullmsg.data(), 1, MPI_INT, 0, world);

        std::vector < char > fullmsg;
        size_t f_size_msg = 0;
        if (!proc_id)
        {
            int i = -1;
            for (auto s : sizefullmsg)
            {
                dispfullmsg[++i] = f_size_msg;
                f_size_msg += static_cast < size_t >( s );
            }
            fullmsg.resize(f_size_msg);
        }
        MPI_Gatherv( msg.data(), size_msg,MPI_CHAR,fullmsg.data(),sizefullmsg.data(), dispfullmsg.data(),MPI_CHAR, 0, world);

        if (!proc_id)
        {
            std::stringstream ss;
            ss<<dim<<"_"<<name;
            std::ofstream out;
            out.open(ss.str().c_str());
            out<<"digraph G {"<<std::endl;
            out.write(fullmsg.data(),f_size_msg);
            out<<"}";
            out.close();

        }
    }

}


} // end namespace
#endif
#endif
