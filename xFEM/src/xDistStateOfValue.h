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

// xlinalg
#include "xDistIndex.h"


namespace xfem
{

//===== xDistStateDofCreator =================================================================================================
/// State dof creator for distributed mesh with non uniform space on it (keys per proc boundary entities may be different).
//! In this class state are only created for key owned by current proc (double manager partition manager gives this information).
//! It respect xStateDofCreator usual conditions (if state exist or test Condition methode is false, no state creation )
//! It respect xStateDofCreator usual numbering (dof use local subset size to creat dof number )
//! The instance of this class must be given to finalizeDof to set a uniform dof numbering across all proc
template < typename Condition = xTrue, typename VT = double >
class xDistStateDofCreator : public Condition
{
public:
    /// create instance
    xDistStateDofCreator(xValueManagerDist<VT> & dm_, const std::string & sub);
    xDistStateDofCreator(const typename xValueManagerDist<VT>::partman_t &part_man_,xValueManagerDist<VT> & dm_, const std::string & sub);
    /// set the state if conditions fulfill
    xValue < VT > * operator()(const xValKey & key, xValue < VT >* v);
    /// give the PM
    const typename xValueManagerDist<VT>::partman_t & getPartitionManager();
    /// give the double manager
    xValueManagerDist<VT> &getValueManager();
    /// give the subset used to store dof
    const std::string &getSubset();
private:
    xValueManagerDist<VT> & dm;
    const typename xValueManagerDist<VT>::partman_t &part_man;
    const std::string subset;
};

//===== xDistStateDofCreatorUniform ==============================================================================================
/// State dof creator for distributed mesh with uniform space on it (keys per proc boundary entities are the same).
//! In this class state are only created for key related to entity owned by current proc.
//! It respect xStateDofCreator usual conditions (if state exist or test Condition method is false, no state creation )
//! It respect xStateDofCreator usual numbering (dof use local subset size to create dof number )
//! The instance of this class must be given to finalizeDofUniform to set a uniform dof numbering across all proc
template < typename PM, typename Condition = xTrue, typename VT = double >
class xDistStateDofCreatorUniform : public Condition
{
public:
    /// create instance
    xDistStateDofCreatorUniform(const PM &part_man_, xValueManagerDist<VT> & val, const std::string & sub);
    /// set the state if conditions fulfill
    xValue < VT > * operator()(const xValKey & key, xValue < VT >* v);
    /// give the PM
    const PM & getPartitionManager();
    /// give the double manager
    xValueManagerDist<VT> &getValueManager();
    /// give the subset used to store dof
    const std::string &getSubset();
private:
    const PM &part_man;
    xValueManagerDist<VT> & val_manager;
    const std::string subset;
};


//===== xDistDofKeyAndInformationManager ============================================================================================
/// Dof information and key manager
//! This class is a key (for information) and information manager used in finalizeDof to exchange dof number.
//
template<typename VT>
class xDistDofKeyAndInformationManager
{
public:
    xDistDofKeyAndInformationManager(xValueManagerDist<VT> & dm_,const typename xValueManagerDist<VT>::partman_t &pm);
    // mandatory types for information class family
    typedef xtool::homogeneous_data_style_trait data_style_trait;
    typedef xtool::send_and_recv_keys_communication_trait communication_trait;
    typedef const xfem::xValKey * information_key_t;
    typedef int information_t;

    // Types related to data that are used to get information key
    // With this class we will iterate on double manager entries to set keys. Data type is then the pair Key/Value
    typedef typename std::iterator_traits < typename xValueManagerDist<VT>::map_iterator >::value_type data_t;

    // mandatory method for information class family
    information_key_t localObjectKey( const data_t &o);
    information_key_t remoteObjectKey(const xtool::xRemoteObject < xValKey > &ro, const data_t &lo );
    xtool::xConstPartitionObject < xValKey > getConstPartitionObject(const data_t &o);
    information_t getInfo(const information_key_t &key, int sendto);
    void setInfo(information_key_t key, const information_t &info, int receivedfrom);

private:
    xValueManagerDist<VT> &dm;
    const typename xValueManagerDist<VT>::partman_t &part_man;
};
//===== xDistDofKeyAndInformationManagerUniform =========================================================================================
/// Dof information and key manager
//! This class is a key (for information) and information manager used in finalizeDofUniform to exchange dof number.
//
template < typename PM, typename VT >
class xDistDofKeyAndInformationManagerUniform
{
public:
    xDistDofKeyAndInformationManagerUniform(const PM &part_man_, xValueManagerDist<VT> & dm_);
    // mandatory types for information class family
    typedef xtool::homogeneous_data_style_trait data_style_trait;
    typedef xtool::send_and_recv_keys_communication_trait communication_trait;
    typedef xfem::xValKey information_key_t;
    typedef int information_t;

    // Types related to data that are used to get information key
    // With this class we will iterate on double manager entries to set keys. Data type is then the pair Key/Value
    typedef typename std::iterator_traits < typename xValueManagerDist<VT>::map_iterator >::value_type data_t;

    // mandatory method for information class family
    information_key_t localObjectKey( const data_t &o);
    information_key_t remoteObjectKey(const xtool::xRemoteObject < AOMD::mEntity > &ro, const data_t &lo );
    xtool::xConstPartitionObject < AOMD::mEntity > getConstPartitionObject(const data_t &o);
    information_t getInfo(const information_key_t &key, int sendto);
    void setInfo(information_key_t key, const information_t &info, int receivedfrom);

private:
    const PM &part_man;
    xValueManagerDist<VT> &dm;
};
//===== xDistFixedKeyManager ============================================================================================
/// Fixed key manager
//! This class is a key (for information)  manager used in finalizeFixed to exchange fixed status to have them uniform across all procs
//

template <typename VT>
class xDistFixedKeyManager
{
public:
    // mandatory traits
    typedef xfem::xValKey information_key_t;

    xDistFixedKeyManager<VT>(xValueManagerDist<VT> & dm_);

    // With this class we will iterate on double manager entries to set keys. Data type is then the pair Key/Value
    typedef typename std::iterator_traits < typename xValueManagerDist<VT>::map_iterator >::value_type data_t;

    // mandatory methods
    information_key_t localObjectKey( const data_t &o);
    std::set < int >  getMessageRanks( const information_key_t &lk);

private:
    xValueManagerDist<VT> &dm;
    const typename xValueManagerDist<VT>::partman_t &part_man;


};
//===== xDistFixedKeyManagerUniform ============================================================================================
/// Fixed key manager
//! This class is a key (for information)  manager used in finalizeFixedUniform to exchange fixed status to have them uniform across all procs
//
template < typename PM, typename VT>
class xDistFixedKeyManagerUniform
{
public:
    // mandatory traits
    typedef xfem::xValKey information_key_t;

    xDistFixedKeyManagerUniform(const PM &part_man_, xValueManagerDist<VT> & dm_);

    // With this class we will iterate on double manager entries to set keys. Data type is then the pair Key/Value
    typedef typename std::iterator_traits < typename xValueManagerDist<VT>::map_iterator >::value_type data_t;

    // mandatory methods
    information_key_t localObjectKey( const data_t &o);
    std::set < int >  getMessageRanks( const information_key_t &lk);

private:
    const PM &part_man;
    xValueManagerDist<VT> &dm;


};
//===== xDistFixedInfoManager ============================================================================================
/// Fixed information manager
//! This class is an information manager used in finalizeFixed to exchange fixed status to have them uniform across all procs
//
template <typename VT>
class xDistFixedInfoManager
{
public:

    // mandatory types
    typedef typename xDistFixedKeyManager<VT>::information_key_t information_key_t;
    typedef std::pair < const xValKey *,VT > information_t;
    typedef xtool::homogeneous_data_style_trait data_style_trait;
    typedef xtool::send_only_keys_communication_trait communication_trait;


    xDistFixedInfoManager( xValueManagerDist<VT> & dm_);
    size_t getNumberOfModification() {return nbc; }

    // mandatory methods
    void setInfo( const std::vector < information_t > &infos, int receivedfrom);
    information_t  getInfo(information_key_t key, int sendto);

private:
    xValueManagerDist<VT> &dm;
    const typename xValueManagerDist<VT>::partman_t &part_man;
    int nbc;
};
//===== xDistFixedInfoManagerUniform ============================================================================================
/// Fixed information manager
//! This class is an information manager used in finalizeFixedUniform to exchange fixed status to have them uniform across all procs
//
template < typename PM, typename VT >
class xDistFixedInfoManagerUniform
{
public:
    typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
    typedef xtool::send_only_keys_communication_trait communication_trait;

    typedef typename xDistFixedKeyManagerUniform < PM, VT >::information_key_t information_key_t;

    xDistFixedInfoManagerUniform(const PM &part_man_, xValueManagerDist<VT> & dm_);
    size_t getNumberOfModification() {return nbc; }

    // mandatory methods
    void getInfo(information_key_t key, xtool::xMpiInputBuffer & buff, int sendto);
    void setInfo(const xtool::xMpiOutputBuffer & buff, int receivedfrom);
    size_t getApproxDataSize() const;

private:
    const PM &part_man;
    xValueManagerDist<VT> &dm;
    int nbc;
};
//===== FinalizeDof ===========================================================================================================
/// FinalizeDof function must be used after all local dof have been created. It finalize state creation and set unique dof id for all procs.
//! For keys which are not already created (key not owned by current proc) a simple xStateDofCreator is used
//! to set there dummy state. A offset is passed to all proc to renumber all local dof. And then a simple scatter exchange correct dof for all procs.
//! This version is only working with double manager partition manager. Says with distributed double manager (all keys nows there counterpart in remotes)
//! This is the most general version compare to FinalizeDofUniform that imply that keys are uniform across process (i.e. on proc boundary entities
//! do have the same set of keys in every process). It implies that double manager partition manager is constructed before.
//
//! Procs involved in this operation are the one corresponding to communicator given by dof_creator argument.
//
//! It returns a shared pointer on a xDistIndex newly allocated instance
//
//
template < typename Condition = xTrue, typename VT >
std::shared_ptr < xlinalg::xDistIndex > finalizeDof( xDistStateDofCreator < Condition, VT > &dof_creator);

//===== FinalizeDofUniform ========================================================================================================
/// FinalizeDofUniform function must be used after all local dof have been created. It finalize state creation and set unique dof id for all procs.
//! For keys which are not already created (key for entity on mesh proc boundary not owned by current proc) a simple xStateDofCreator is used
//! to set there dummy state. A offset is passed to all proc to renumber all local dof. And then a simple scatter exchange correct dof for all procs
//! This is less general than FinalizeDof. It consider that keys are uniform across process (i.e. on proc boundary entities
//! do have the same set of keys in every process). Compare to FinalizeDof, double manager partition manager need not to be created. The benefit is then
//! the gain in communication/computation time by not creating this partition.
//
//! Procs involved in this operation are the one corresponding to communicator given by dof_creator argument.
//
//! It returns a shared pointer on a xDistIndex newly allocated instance
//
template < typename PM, typename Condition = xTrue, typename VT >
std::shared_ptr < xlinalg::xDistIndex > finalizeDofUniform( xDistStateDofCreatorUniform < PM, Condition, VT > &dof_creator);

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
//! Compare to FinalizeFixedUniform this version act directely with help of double manager partition manager that must have
//! been created by genPartitionManager. If it is not the case an assert will stop this function when internaly it will try
//! to retrive partition manager from dm.
//! By using dm partition manager this function will be able to deal with the case of key existing only in one proc and not in
//! others for some proc boundary entities.
//
//! Procs involved in this operation are the one corresponding to communicator of double manager partition manager.
//
//! It return a bool which is true only if this function change local status to fixed
//
template<typename VT>
bool finalizeFixed( xValueManagerDist<VT> & dm){
    // create key manager
    xfem::xDistFixedKeyManager<VT> key_info_manager(dm);

    // set keys for information to be exchanged
    xtool::xKeyContainerSendOrRecv < typename xfem::xDistFixedKeyManager<VT>::information_key_t, xfem::xHashValKey, xfem::xEqualValKey > infokey_container(dm.getPartitionManager().getComm());
    infokey_container.accumulateKeys(dm.begin(),dm.end(),key_info_manager);

    // create info manager
    xfem::xDistFixedInfoManager<VT> info_manager(dm);

    // exchange information (create fixed status if needed)
    exchangeInformation(infokey_container,info_manager);

    return ( info_manager.getNumberOfModification() != 0 );
}

//===== FinalizeFixedUniform =========================================================================================================
/// FinalizeFixedUniform function must be used after all local fixed values have been set. It finalize state creation among all procs.
//! Compare to FinalizeFixed this version works only if keys are uniform across process (i.e. on proc boundary, entities
//! do have the same set of keys in every process).
//! Compare to finalizeFixed, double manager partition manager need not to be created (may be a benefit)
//
//
//! Procs involved in this operation are the one corresponding to communicator of double manager given as argument.
//
//! It return a bool which is true only if this function change local status to fixed
//
template < typename PM, typename VT >
bool finalizeFixedUniform(  const PM &part_man,  xValueManagerDist<VT> & dm);

}  // end namespace

#include "xDistStateOfValue_imp.h"








#endif
