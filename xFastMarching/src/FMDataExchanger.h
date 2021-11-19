/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */
#ifndef _FMDATAEXCHANGER_
#define _FMDATAEXCHANGER_


#include "xDataExchanger.h"
#include "xMPIDataType.h"


// ==================================================================================
// Complete  Transport concept requirements when in parallel
// Required function to pack/unpack transport concept object
// must be out of namespace to avoid any namespace artificial dependency. It  stay general
template < class SCAL >
inline void pack(const xfastmarching::transportConcept < SCAL > &trans,xtool::xMpiInputBuffer &buff) {; }
template < class SCAL >
inline void unPack(xfastmarching::transportConcept < SCAL > &trans, const xtool::xMpiOutputBuffer &buff){; }
// ==================================================================================

namespace xfastmarching
{
namespace internal
{


template < typename MI, typename UPDATER >
class FMDataExchanger
{
    public:
        typedef typename MI::vertex vertex;
        typedef typename MI::entity entity;
        typedef typename MI::partman_t partman_t;
        typedef typename MI::key_manager_sor_t key_manager_sor_t;
        typedef typename UPDATER::scal scal;
        // mandatory traits for info manager
        typedef typename MI::information_key_t information_key_t;
        typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
        typedef xtool::send_only_keys_communication_trait communication_trait;
        typedef std::pair < scal, const vertex * > pair_t;
        typedef std::unordered_map < const vertex *, scal > trial_t;

        FMDataExchanger(const MI & mi_,UPDATER &updater_, trial_t & bTrial_) :
            mi(mi_)
            ,updater(updater_)
            ,bTrial(bTrial_)
            ,pm(mi.reg.getPartitionManager())
            ,key_man(pm)
            ,key_container(pm.getComm())
            ,L(updater.L)
            ,min_exchanged(L)
        {
            setFilteredRegion(mi);
            updater.initBoundary(beginBnd(mi),endBnd(mi));

        }

        inline MPI_Comm  getComm() const { return pm.getComm(); }
        inline void  clearKeyContainer() { key_container.clearKeys(); }

        scal exchange()
        {
            min_exchanged = L;
            xtool::exchangeInformation(key_container,*this);
            return min_exchanged;
        }
        inline void accumulate(const vertex &vk)
        {
            if ( updater.isBoundaryLocalyUpdated (vk) )
                key_container.accumulateKeys(vk,key_man);
        }

        // mandatory methods for communication
        void setInfo(const xtool::xMpiOutputBuffer & buff, int receivedfrom)
        {

            while (!buff.exhausted() )
            {
                auto trial_pair = updater.unPack(buff,receivedfrom);
                const vertex *v = trial_pair.second;
                if (v)
                {
                    /*
                       scal Tv = trial_pair.first;
                       if (Tv < min_exchanged) min_exchanged = Tv;
                       bTrial[v]=Tv;
                     */
                    bTrial[v] = trial_pair.first;
                    scal Tv;
                    updater.getGlobalLs(*v,Tv);
                    if (Tv < min_exchanged) min_exchanged = Tv;
                }
            }
        }
        void getInfo(information_key_t key, xtool::xMpiInputBuffer & buff, int sendto)
        {
            if (updater.status(static_cast < const vertex * >( key )))
            {
                const vertex *v = static_cast < const vertex * >( key );

                const vertex*  vr = getRemoteVertex(mi, *v,sendto);
                buff.pack(&vr,1,MPI_AINT);
                updater.pack(*v,buff);
            }
        }
        size_t getApproxDataSize(void)
        {
            return sizeof( scal )+sizeof( vertex * );
        }

    private:
        const MI &mi;
        UPDATER &updater;
        trial_t & bTrial;
        const partman_t &pm;
        key_manager_sor_t key_man;
        xtool::xKeyContainerSendOrRecv < information_key_t >  key_container;
        const scal L;
        scal min_exchanged;
};

} // namespace
} // namespace

#endif
