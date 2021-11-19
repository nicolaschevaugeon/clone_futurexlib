/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "xValueManager.h"
#include "xDataExchanger.h"

namespace xfem
{

// =======================================================================================================================================
// Local class to exchange information
class xDistKeyKeymanager
{
    public:
        // mandatory traits
        typedef xValKey information_key_t;

        xDistKeyKeymanager(const xMesh::partman_t &part_man_) :
            part_man(part_man_)
        {}

        // With this class we will iterate on double manager entries to set keys. Data type is then the pair Key/Value
        typedef std::iterator_traits < xValueManagerDist<double>::map_iterator >::value_type data_t;

        // mandatory methods
        information_key_t localObjectKey( const data_t &o)
        {
            return o.first;
        }
        std::set < int >  getMessageRanks( const information_key_t &key)
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
                // grabe remote proc id for object and store them in res
                for (const auto ro : po.getRemoteObjectsCollectionRange( ))
                {
                    res.insert(ro.getProcId());
                }
            }

            return res;

        }

    private:
        const xMesh::partman_t &part_man;
};
// this structure is here to transfert the xvalkey information with its origine
struct transXvalKey
{
    xValKey * adress;
    AOMD::mEntity*  recever_Enti;
    int Refe;
    short Phys;
    short Geom;
};
class xDistKeyInfoManager
{
    public:
        typedef xtool::homogeneous_data_style_trait data_style_trait;
        typedef xtool::send_only_keys_communication_trait communication_trait;


        typedef xDistKeyKeymanager::information_key_t information_key_t;
        typedef transXvalKey information_t;


        xDistKeyInfoManager(xValManagerKeyPartitionManager < xValKey,  xValue < double >, xHashValKey, xEqualValKey, xMesh::partman_t >::partman_t &part_man_, const xValueManagerDist<double> & dm_,const xMesh::partman_t &mesh_partman_) :
            part_man(part_man_)
            ,dm(dm_)
            ,mesh_part_man(mesh_partman_)
        {}

        // mandatory methods
        information_t getInfo(const information_key_t &key, int sendto)
        {
            // Get adresse in dm of the key
            const xValKey *lka=dm.findKeyAdress(key);
            assert(lka);

            // first get entity of this key to create remote
            AOMD::mEntity *lo = key.getEnti();

            // get remote object
            xtool::xConstPartitionObject < AOMD::mEntity > po = mesh_part_man.getConstPartitionObject(*lo);

            assert (po.hasRemoteObject());
            transXvalKey tmp;
            tmp.adress = const_cast<xValKey *>(lka);
            tmp.recever_Enti = const_cast<AOMD::mEntity *>(po.getRemoteObjectOn(sendto));
            tmp.Refe = key.getRefe();
            tmp.Phys = key.getPhys();
            tmp.Geom = key.getGeom();
            return tmp;
        }
        void setInfo(const std::vector < information_t > &infos, int receivedfrom )
        {
            for (auto info : infos)
            {
                xValKey tmp(info.Phys,info.Geom,info.recever_Enti,info.Refe);
                
                // if this key is present in dm then a connection must be set
                if (dm.findKeyAdress(tmp))
                {
                    auto npo = part_man.getPartitionObject(tmp);
                    npo.insert(receivedfrom ,info.adress);
                }
            }
        }
    private:
        xValManagerKeyPartitionManager < xValKey,  xValue < double >, xHashValKey, xEqualValKey, xMesh::partman_t >::partman_t &part_man;
        const xValueManagerDist<double> &dm;
        const xMesh::partman_t &mesh_part_man;
};



// =======================================================================================================================================
// xValManagerKeyPartitionManager specialization implementation

xValManagerKeyPartitionManager < xValKey,xValue < double >,xHashValKey,xEqualValKey,xMesh::partman_t >::xValManagerKeyPartitionManager(const xValueManagerDist<double> &dm_,const xMesh::partman_t & mesh_partman_) :
    dm(dm_)
    ,mesh_partman(mesh_partman_)
    ,partition_manager(mesh_partman.getComm())
{}

/// function that create distributed informations
void xValManagerKeyPartitionManager < xValKey,xValue < double >,xHashValKey,xEqualValKey,xMesh::partman_t >::genPartitionManager(void)
{
    partition_manager.clear();
    MPI_Comm world = partition_manager.getComm();


    // set keys: xValkey that are on some proc boundary
    xDistKeyKeymanager key_manager(mesh_partman);
    xtool::xKeyContainerSendOrRecv < xDistKeyKeymanager::information_key_t, xHashValKey, xEqualValKey > key_container(world);
    key_container.accumulateKeys(dm.begin(),dm.end(),key_manager);

    // exchange: creation of partition_manager
    xDistKeyInfoManager info_manager(partition_manager,dm,mesh_partman);
    xtool::exchangeInformation(key_container,info_manager);

}
auto xValManagerKeyPartitionManager < xValKey,xValue < double >,xHashValKey,xEqualValKey,xMesh::partman_t >::getPartitionManager(void) const->const partman_t &
{
    return partition_manager;
}

} // end of namespace
