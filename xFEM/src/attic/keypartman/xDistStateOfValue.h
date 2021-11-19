/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef XDISTSTATEOFVALUE_H
#define XDISTSTATEOFVALUE_H

//std
#include <cassert>
#include <string>


// parrallel
#include <mpi.h>

// xfem
#include "xValueManager.h"
#include "xStateOfValue.h"

// eXlibris_tools
#include "xPartitionManager.h"
#include "xDataExchanger.h"

// lalg
#include "xDistIndex.h"


namespace xfem
{

//===== xDistStateDofCreator =================================================================================================
/// State dof creator for distributed mesh
//! In this class state are only created for key related to entity owned by curent proc.
//! It respect xStateDofCreator usual conditions (if state exist or test Condition methode is false, no state creation )
//! It respect xStateDofCreator usual numbering (dof use local subset size to creat dof number )
//! The instance of this class must be given to finalizeDof to set a uniform dof numbering across all proc
template < typename PM, typename Condition = xTrue >
class xDistStateDofCreator : public Condition
{
    public:
        /// create instance
        xDistStateDofCreator(const PM &part_man_, xValueManagerDist<double> & val, const std::string & sub);
        /// set the state if conditions fulfill
        xValue < double > * operator()(const xValKey & key, xValue < double >* v);
        /// give the PM
        const PM & getPartitionManager();
        /// give the double manager
        xValueManagerDist<double> &getDoubleManager();
        /// give the subset used to store dof
        const std::string &getSubset();
    private:
        const PM &part_man;
        xValueManagerDist<double> & val_manager;
        const std::string subset;
};

//===== xDistStateDofCreatorOLD ==============================================================================================
/// State dof creator for distributed mesh
//! In this class state are only created for key related to entity owned by curent proc.
//! It respect xStateDofCreator usual conditions (if state exist or test Condition methode is false, no state creation )
//! It respect xStateDofCreator usual numbering (dof use local subset size to creat dof number )
//! The instance of this class must be given to finalizeDof to set a uniform dof numbering across all proc
template < typename PM, typename Condition = xTrue >
class  xDistStateDofCreatorOLD : public Condition
{
    public:
        /// create instance
        xDistStateDofCreatorOLD(const PM &part_man_, xValueManagerDist<double> & val, const std::string & sub);
        /// set the state if conditions fulfill
        xValue < double > * operator()(const xValKey & key, xValue < double >* v);
        /// give the PM
        const PM & getPartitionManager();
        /// give the double manager
        xValueManagerDist<double> &getDoubleManager();
        /// give the subset used to store dof
        const std::string &getSubset();
    private:
        const PM &part_man;
        xValueManagerDist<double> & val_manager;
        const std::string subset;
};


//===== xDistDofKeyAndInformationManager ============================================================================================
/// Dof information and key manager
//! This class is a key (for information) and information manager used in finalizeDof to exchange dof number.
//
template < typename PM >
class xDistDofKeyAndInformationManager
{
    public:
        xDistDofKeyAndInformationManager(const PM &part_man_, xValueManagerDist<double> & dm_);
        // mandatory types for information class family
        typedef xtool::homogeneous_data_style_trait data_style_trait;
        typedef xtool::send_and_recv_keys_communication_trait communication_trait;
        typedef const xfem::xValKey * information_key_t;
        typedef int information_t;

        // Types related to data that are used to get information key
        // With this class we will iterate on double manager entries to set keys. Data type is then the pair Key/Value
        typedef std::iterator_traits < xValueManagerDist<double>::map_iterator >::value_type data_t;

        // mandatory method for information class family
        information_key_t localObjectKey( const data_t &o);
        information_key_t remoteObjectKey(const xtool::xRemoteObject < xValKey > &ro, const data_t &lo );
        xtool::xConstPartitionObject < xValKey > getConstPartitionObject(const data_t &o);
        information_t getInfo(const information_key_t &key, int sendto);
        void setInfo(information_key_t key, const information_t &info, int receivedfrom);

    private:
        const PM &part_man;
        xValueManagerDist<double> &dm;
};
//===== xDistDofKeyAndInformationManagerOLD =========================================================================================
/// Dof information and key manager
//! This class is a key (for information) and information manager used in finalizeDof to exchange dof number.
//
template < typename PM >
class xDistDofKeyAndInformationManagerOLD
{
    public:
        xDistDofKeyAndInformationManagerOLD(const PM &part_man_, xValueManagerDist<double> & dm_);
        // mandatory types for information class family
        typedef xtool::homogeneous_data_style_trait data_style_trait;
        typedef xtool::send_and_recv_keys_communication_trait communication_trait;
        typedef xfem::xValKey information_key_t;
        typedef int information_t;

        // Types related to data that are used to get information key
        // With this class we will iterate on double manager entries to set keys. Data type is then the pair Key/Value
        typedef std::iterator_traits < xValueManagerDist<double>::map_iterator >::value_type data_t;

        // mandatory method for information class family
        information_key_t localObjectKey( const data_t &o);
        information_key_t remoteObjectKey(const xtool::xRemoteObject < AOMD::mEntity > &ro, const data_t &lo );
        xtool::xConstPartitionObject < AOMD::mEntity > getConstPartitionObject(const data_t &o);
        information_t getInfo(const information_key_t &key, int sendto);
        void setInfo(information_key_t key, const information_t &info, int receivedfrom);

    private:
        const PM &part_man;
        xValueManagerDist<double> &dm;
};
//===== xDistFixedKeyManager ============================================================================================
/// Fixed key manager
//! This class is a key (for information)  manager used in finalizeFixed to exchange fixed status to have them uniform across all procs
//
template < typename PM >
class xDistFixedKeyManager
{
    public:
        // mandatory traits
        typedef xfem::xValKey information_key_t;

        xDistFixedKeyManager(const PM &part_man_, xValueManagerDist<double> & dm_);

        // With this class we will iterate on double manager entries to set keys. Data type is then the pair Key/Value
        typedef std::iterator_traits < xValueManagerDist<double>::map_iterator >::value_type data_t;

        // mandatory methods
        information_key_t localObjectKey( const data_t &o);
        std::set < int >  getMessageRanks( const information_key_t &lk);

    private :
        const PM &part_man;
        xValueManagerDist<double> &dm;


};
//===== xDistFixedInfoManager ============================================================================================
/// Fixed information manager
//! This class is an information manager used in finalizeFixed to exchange fixed status to have them uniform across all procs
//
template < typename PM >
class xDistFixedInfoManager
{
    public:
        typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
        typedef xtool::send_only_keys_communication_trait communication_trait;

        typedef typename xDistFixedKeyManager < PM >::information_key_t information_key_t;

        xDistFixedInfoManager(const PM &part_man_, xValueManagerDist<double> & dm_);
        size_t getNumberOfModification() {return nbc;}

        // mandatory methods
        void getInfo(information_key_t key, xtool::xMpiInputBuffer & buff, int sendto);
        void setInfo(const xtool::xMpiOutputBuffer & buff, int receivedfrom);
        size_t getApproxDataSize() const;

    private :
        int nbc;
        const PM &part_man;
        xValueManagerDist<double> &dm;
};
//===== FinalizeDof ===========================================================================================================
/// FinalizeDof function must be used after all local dof have been created. It finalize state creation and set unique dof id for all procs.
//! For keys which are not already created (key for entity on mesh proc boundary not owned by current proc) a simple xStateDofCreator is used
//! to set there dummy state. A offset is passed to all proc to renumber all local dof. And then a simple scatter exchange correct dof for all procs
//
//! Procs involved in this operation are the one corresponding to communicator given by dof_creator argument.
//
//! It returns a shared pointer on a xDistIndex newly allocated instance
//
template < typename PM, typename Condition = xTrue >
std::shared_ptr < lalg::xDistIndex > finalizeDof( xDistStateDofCreator < PM, Condition > &dof_creator);

//===== FinalizeDofOLD ========================================================================================================
/// FinalizeDof function must be used after all local dof have been created. It finalize state creation and set unique dof id for all procs.
//! For keys which are not already created (key for entity on mesh proc boundary not owned by current proc) a simple xStateDofCreator is used
//! to set there dummy state. A offset is passed to all proc to renumber all local dof. And then a simple scatter exchange correct dof for all procs
//
//! Procs involved in this operation are the one corresponding to communicator given by dof_creator argument.
//
//! It returns a shared pointer on a xDistIndex newly allocated instance
//
template < typename PM, typename Condition = xTrue >
std::shared_ptr < lalg::xDistIndex > finalizeDofOLD( xDistStateDofCreatorOLD < PM, Condition > &dof_creator);

//===== FinalizeFixed =========================================================================================================
/// FinalizeFixed function must be used after all local fixed values have been set. It finalize state creation among all procs.
//! Here is a example in 3D which explain why we need to use this function.
//! Say user fixed values via a part surface condition. A proc may have region element connect to this surface by different
//! ways:
//!         It may be connect by face of its element which will naturally be treated by TreatmentOfEssentialEnv family
//!         It may be connect by edge or node of its element which will not be treated by TreatmentOfEssentialEnv family
//! If two adjacent proc domain have respectively all there element in both category, local double manager will have different
//! status for value common to both domain and on the imposed surface. For domain with face on surface, common value will be
//! fixed. And for domain with edge or node on surface  common value will be dof. It can't be !
//! Fixed status must be consistent across proc...
//
//! Procs involved in this operation are the one corresponding to communicator given by argument.
//
//! It return a bool which is true only if this function change local status to fixed
//
template < typename PM >
bool finalizeFixed(  const PM &part_man,  xValueManagerDist<double> & dm);

}  // end namespace

#include "xDistStateOfValue_imp.h"

#endif
