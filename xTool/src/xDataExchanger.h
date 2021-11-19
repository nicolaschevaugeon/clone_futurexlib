/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef XDATAEXCHANGER_H
#define XDATAEXCHANGER_H

//std
#include <cassert>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <unordered_set>
#include <string>
#include <exception>
#include <sstream>
#include <type_traits>

// parrallel
#include <mpi.h>

//eXlibris_tools
#include "xDataExchangerTraits.h"
#include "xPartitionManager.h"
#include "xMPITag.h"

/*! \file
 *  \brief Function that help exchanging data in between process
 *
 * General concept
 * ---------------
 * This file provide function exchangeInformation which use the concept of 'associated information' to
 * transfer information in between process.
 * Associated information imply that a key (unique among process) exist for any information that we want
 * to transfer. This set of keys is used to retrieve associated information (that can be of any kind).
 * Those information are pack and/or unpack in/from communication buffer before/after communication.
 *
 * Compare to a simple buffer that is send or receive with MPI function this function offers generality
 * and performance :
 * <ul><li>  When key exist at both end of the communicating process, sorting keys avoid
 * transferring keys themselves. This is efficient especially when information is smaller then the key.
 * </li><li> A set of key may be create once and used many times for many type of informations.
 * </li><li> Heterogeneous or homogeneous information may be used.
 * </li><li> When key are only known in sending process, implementation use MPI3 norm to exchange
 * information in a rather efficient way. It provide a scalable solution to the
 * "Dynamic Sparse Data Exchange" problem (see Scalable communication protocols for dynamic sparse
 * data exchange, T. Hoefler et al., Proceedings of the 15th ACM SIGPLAN Symposium on Principles
 * and Practice of Parallel Programming Pages 159-168 )
 * </li><li> In most other communication pattern, posting non blocking receive are used to overlap
 * communication and computation (at least unpacking computation)
 * </li><li> A set of key may be create once and used indifferently to send or receive information
 * from other proc (pattern when you first send new information from a source to a target and then
 * receive in the source a "answer" from the target ).
 * </li><li> The separation of the information from the key that reference it permit orthogonalization
 * of implementation and, hopefully, permit future extensions (to windowing communication for example)
 * </li><li> Packing and unpacking thanks to the use of dedicated manager (information manager) may be
 * inserted in the heart of algorithm. Says no intermediate buffers are expected and computation are likely
 * to be done during those operations.
 * </li></ul>
 *
 * Template parameter name dictionary
 * ----------------------------------
 *
 * <ul><li>    K    refer to the type of key type
 * </li><li>   KC   refer to the type of key container type
 * </li><li>   KM   refer to the type of key manager type
 * </li><li>   IM   refer to the type of information manager type
 * </li><li>   ITER refer to the type of iterator on object or on pointer to object
 * </li><li>   T    refer to the type of object given by ITER iterator
 * </li><li>   O    refer to the type of object used in xPartitionManager. It discribe element (node,edge,face,index, adress, ...) that are 
 *                  duplicate in an other process and managed by a xPartitionManager. O may be T or not. For example T may be an xValKey
 *                  and O a mEntity (see xFem library)
 * </li><li>   I    refer to the type of information type 
 * </li></ul>
 *
 *
 * Concept of a key container
 * --------------------------
 * The keys have to be stored in appropriate containers (actually class xKeyContainerSendAndRecv and
 * xKeyContainerSendOrRecv are providing rather generic storage fulfilling requirement of
 * exchangeInformation ).
 *
 * Those key information container must provide:
 * <ul><li>  information_key_t type: it describe the nature of the keys that are used to access information
 * </li><li> keys_association_trait type: it is used to select correct exchanging policy considering that key
 * are know on both communication side (send_and_recv_keys_association_trait) or only know on one
 * communication side (send_or_recv_keys_association_trait, on the sender side or on the receiver side)
 * </li><li> acces (via friend token) to  class exchangeInformationPolicy
 * </li><li> for send_and_recv_keys_association_trait:
 * <ul><li> A STL sorted associative container send_keys_to that store for a proc id to communicate with, its set of
 * pair of sorted key (the key in local proc and the corresponding key in remote process).
 * </li><li> A STL sorted associative container recv_keys_from that store for a proc id to communicate with, its set of sorted key.
 * </li></ul>
 * </li><li> for send_or_recv_keys_association_trait:
 * <ul><li> A STL sorted associative container keys_from_to that store for a proc id to communicate with, its set of unsorted key.
 * </li></ul>
 * </li><li> accumulate method to store keys using iterator (ITER type) and key manager (KM type)
 * </li><li> clearKeys method to clear containers.
 * </li><li> getComm method to give communicator related to container.
 * </li></ul>
 *
 * The idea is to iterate on a collection of object (o) that are related to information that have to be transfered (via accumulate method). The key manager (usually
 * through a xPArtitionManager member) is use to choose to create (or not)  an object o key  which is stored then in the key container.
 *
 * Accumulate method provided by xKeyContainerSendAndRecv are:
 * <ul><li>  accumulateKeysOwnerScatter: message for selected key (by the key manager) are expected to be send
 * from owner (see xPartitionManager.h) to all remote process.
 * </li><li> accumulateKeysOwnerGather: message for selected key (by the key manager) are expected to be retrieved
 * on owner from all remote process.
 * </li><li> accumulateKeysAllGather: message for selected key (by the key manager) are expected to be send and receive
 * to all remote process and from all remote process.
 * </li></ul>
 *
 *
 * Accumulate method provided by xKeyContainerSendOrRecv is accumulateKeys: message for selected key (by the key manager)
 * are  expected to be send or receive to remote process given by user via key manager.
 *
 * Accumulate method may be used many time. The current state of the containers is taken into account by exchangeInformation.
 *
 * Concept of key manager
 * ----------------------
 *
 * The key manager used by accumulate method (of key container class) must provide :
 * <ul><li>  for xKeyContainerSendAndRecv :
 * <ul><li>  A localObjectKey method that take an instance to the local object (lo) investigate as argument. It must return
 * local key, the key associated to this local object.
 * </li><li>  A getConstPartitionObject  method that take lo as argument and return a ConstPartitionObject po (see xPartitionManager.h)
 * </li><li>  A remoteObjectKey method that take ro and lo as argument. ro is the remote object of lo given by po. It return
 * a remote key rk that is used for sorting key set in order of remote process.
 * </li></ul>
 * </li><li>  for xKeyContainerSendOrRecv :
 * <ul><li>  A localObjectKey method that take an instance to the local object (lo) investigate as argument. It must return
 * local key (lk), the key associated to this local object.
 * </li><li>  A getMessageRanks  method that take lk as argument and return a set of int, giving process rank to whom communication
 * is expected
 * </li></ul>
 * </li><li>  information_key_t type: it describe the nature of the keys that are used to access information. It must be equal to
 * the one provided by class calling accumulate method.
 * </li></ul>
 *
 * User is responsible of key manger creation. It is given to accumulate method so that correct choice (or filtering) is done
 * to select keys. Those keys are expected to be associated to information to communicate.
 *
 * Key manager class required structure summary
 * --------------------------------------------
 *
 * For xKeyContainerSendAndRecv :
 * \code{.cpp}
 * template < typename K,  typename O, typename T >
 * class KeyManagerConcept
 * {
 *  public:
 *     // mandatory traits
 *     typedef K information_key_t;
 *
 *     // mandatory methods
 *     information_key_t localObjectKey( const T & o);
 *     information_key_t remoteObjectKey(const xtool::xRemoteObject < O > & ro, const T & lo);
 *     xtool::xConstPartitionObject < O > getConstPartitionObject( const T & lo);
 * };
 * \endcode
 *
 * For xKeyContainerSendOrRecv :
 * \code{.cpp}
 * template < typename K,  typename T >
 * class KeyManagerConcept
 * {
 *  public:
 *     // mandatory traits
 *     typedef K information_key_t;
 *
 *     // mandatory methods
 *     information_key_t localObjectKey( const T & o);
 *     std::set < int > getMessageRanks( const information_key_t & lk);
 * };
 * \endcode
 *
 * Concept of info manager
 * -----------------------
 *
 *  This is a user feature that drives information management in exchangeInformation function.
 *  It must provide :
 * <ul><li>   A communication_trait type. It control exchange nature :
 * <ul><li>  send_and_recv_keys_communication_trait: information are sent and receive among all process involved. It must be
 * used only with key container of type send_and_recv_keys_association_trait (xKeyContainerSendAndRecv for now)
 * </li><li>  send_only_keys_communication_trait:  information managed by keys are sent to involved process. It must be
 * used only with key container of type send_or_recv_keys_association_trait (xKeyContainerSendOrRecv for now)
 * </li><li>  recv_only_keys_communication_trait:  information managed by keys are received by involved process. It must be
 * used only with key container of type send_or_recv_keys_association_trait (xKeyContainerSendOrRecv for now)
 * </li></ul>
 * </li><li>  A data_style_trait type. It control information nature :
 * <ul><li>  homogeneous_data_style_trait: information are identical in type.
 * </li><li>  nonhomogeneous_data_style_trait:  information are not identical in type. User is using xMpiInputBuffer and
 * xMpiOutputBuffer class to pack/unpack information of different size.
 * </li></ul>
 * </li><li>  A information_key_t type. It must be equal to the one provided by key_info_container argument of  exchangeInformation function.
 * It correspond to the type used for key.
 * </li></ul>
 *
 * Depending on those trait the following must be provided
 * <ul><li>   send_and_recv_keys_communication_trait/homogeneous_data_style_trait:
 * <ul><li>   A getInfo method to retrieve information associated to the key given as argument for target process
 * given also as argument. The return value is the information which is copy in communication buffer of exchangeInformation policy.
 * </li><li>  A  setInfo method to set  information associated to  the key given as argument for target process
 * given also as argument. The information to set is provided as argument to this function.
 * </li><li>  A  information_t type corresponding to information type
 * </li></ul>
 * </li><li>  send_and_recv_keys_communication_trait/nonhomogeneous_data_style_trait
 * <ul><li>   A getInfo method to retrieve information associated to the key given as argument for target process
 * given also as argument. The information is packed in communication buffer of exchangeInformation policy by using a extra argument of type
 * xMpiInputBuffer. This instance may be used the way user want in this function and it is his/her responsability to provide a way to unpack
 * correctly data in setInfo method.
 * </li><li>  A  setInfo method to set  information associated to the key given as argument for target process
 * given also as argument. The information to set is provided as an xMpiOutputBuffer argument type to this function. Has said above it is user
 * responsibility to unpack correctly information related to keys in this function.
 * </li><li>  A  getApproxDataSize method to estimate mean size in byte (sizeof(char)) of individual information associated to
 * a key that will be send.
 * </li></ul>
 * </li><li>  send_only_keys_communication_trait/homogeneous_data_style_trait
 * <ul><li>   A getInfo method to retrieve information associated to the key given as argument for target process
 * given also as argument. The return value is the information which is copy in communication buffer of exchangeInformation policy.
 * </li><li>  A  setInfo method to set  information received from source process given as argument. The information to set is provided as
 * a STL vector (of information_t values) argument to this function.
 * </li><li>  A  information_t type corresponding to information type
 * </li></ul>
 * </li><li>  send_only_keys_communication_trait/nonhomogeneous_data_style_trait
 * <ul><li>   A getInfo method to retrieve information associated to the key given as argument for target process
 * given also as argument. The information is packed in communication buffer of exchangeInformation policy by using a extra argument of type
 * xMpiInputBuffer. This instance may be used the way user want in this function and it is his/her responsibility to provide a way to unpack
 * correctly data in setInfo method.
 * </li><li>  A  setInfo method to set  information  received from source process given as argument. The information to set is provided as
 * an xMpiOutputBuffer argument type to this function. Has said above it is user responsibility to unpack correctly information.
 * </li><li>  A  getApproxDataSize method to estimate mean size in byte (sizeof(char)) of individual information associated to
 * a key that will be send.
 * </li></ul>
 * </li><li>  recv_only_keys_communication_trait/homogeneous_data_style_trait
 * <ul><li>   A getInfo method to retrieve information for target process  given as argument. The return value is a STL vector (of information_t values).
 * </li><li>  A  setInfo method to set  information associated to  the key given as argument for target process
 * given also as argument. The information to set is provided as argument to this function.
 * </li><li>  A  information_t type corresponding to information type
 * </li></ul>
 * </li><li>  recv_only_keys_communication_trait/nonhomogeneous_data_style_trait: not (yet) implemented
 * </li></ul>
 *
 *
 * Information manager class required structure summary
 * ----------------------------------------------------
 *
 * For send_and_recv_keys_communication_trait/homogeneous_data_style_trait :
 * \code{.cpp}
 * template < typename K, typename I >
 * class InfoManagerConcept
 * {
 *  public:
 *     // mandatory traits
 *     typedef K information_key_t;
 *     typedef I information_t;
 *     typedef xtool::homogeneous_data_style_trait data_style_trait;
 *     typedef xtool::send_and_recv_keys_communication_trait communication_trait;
 *
 *     // mandatory methods
 *     void setInfo(information_key_t key, const information_t &info, int receivedfrom);
 *     information_t getInfo(information_key_t key, int sendto);
 * };
 * \endcode
 *
 * For send_and_recv_keys_communication_trait/nonhomogeneous_data_style_trait :
 * \code{.cpp}
 * template < typename K  >
 * class InfoManagerConcept
 * {
 *  public:
 *     // mandatory traits
 *     typedef K information_key_t;
 *     typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
 *     typedef xtool::send_and_recv_keys_communication_trait communication_trait;
 *
 *     // mandatory methods
 *      void setInfo(information_key_t key, const xtool::xMpiOutputBuffer & buff , int receivedfrom);
 *      void getInfo(information_key_t key, xtool::xMpiInputBuffer & buff, int sendto);
 *      size_t getApproxDataSize(void);
 * };
 * \endcode
 *
 * For send_only_keys_communication_trait/homogeneous_data_style_trait :
 * \code{.cpp}
 * template < typename K, typename I >
 * class InfoManagerConcept
 * {
 *  public:
 *     // mandatory traits
 *     typedef K information_key_t;
 *     typedef I information_t;
 *     typedef xtool::homogeneous_data_style_trait data_style_trait;
 *     typedef xtool::send_only_keys_communication_trait communication_trait;
 *
 *     // mandatory methods
 *     void setInfo( const std::vector < information_t > &infos, int receivedfrom);
 *     information_t getInfo(information_key_t key, int sendto);
 * };
 * \endcode
 *
 * For send_only_keys_communication_trait/nonhomogeneous_data_style_trait :
 * \code{.cpp}
 * template < typename K >
 * class InfoManagerConcept
 * {
 *  public:
 *     // mandatory traits
 *     typedef K information_key_t;
 *     typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
 *     typedef xtool::send_only_keys_communication_trait communication_trait;
 *
 *     // mandatory methods
 *      void setInfo(const xtool::xMpiOutputBuffer & buff , int receivedfrom);
 *      void getInfo(information_key_t key, xtool::xMpiInputBuffer & buff, int sendto);
 *      size_t getApproxDataSize(void);
 * };
 * \endcode
 *
 * For recv_only_keys_communication_trait/homogeneous_data_style_trait :
 * \code{.cpp}
 * template < typename K, typename I >
 * class InfoManagerConcept
 * {
 *  public:
 *     // mandatory traits
 *     typedef K information_key_t;
 *     typedef I information_t;
 *     typedef xtool::homogeneous_data_style_trait data_style_trait;
 *     typedef xtool::recv_only_keys_communication_trait communication_trait;
 *
 *     // mandatory methods
 *     void setInfo(information_key_t &key, const information_t &info, int receivedfrom);
 *     std::vector < information_t > getInfo( int sendto);
 * };
 * \endcode
 *
 */

namespace xtool
{


/*!
 * \brief Function which do an exchange of information for keys given by key_info_container following protocol given by info_manager
 * \param [in] key_info_container object that store keys that give a way to retrieve information to exchange. This instance type is one of the two
 *                                class xKeyContainerSendAndRecv or xKeyContainerSendOrRecv. Or it could be any class you want that follows concept
 *                                of key container.
 * \param [in] info_manager object that manage information associated to keys of key_info_container.
 * \param [in] work a function that may be provided to overlap computation and communication. In exchangeInformation when communication is done call
 *  to work function may be done once or many time during incoming message pooling.
 */
template < typename KC, typename  IM >
void exchangeInformation( KC & key_info_container, IM & info_manager,  void ( *work )() = nullptr);




/// Container  that store the Keys in a send and receive keys association manner. Template on K, the type of key.
template < typename K >
class xKeyContainerSendAndRecv
{
    public:
        typedef  K information_key_t;
        typedef send_and_recv_keys_association_trait keys_association_trait;

        xKeyContainerSendAndRecv( MPI_Comm _mpi_comm);

        /// each local object in the iterator range send to all it's remote object if the object is owned by the current proc. the data send to the remote object is not necesseraly the same.
        /// Collective on the MPI_Comm of the exchanger
        /// Only avalaible for exchangeInformationPolicy with trait send_recv_keys_association_trait
        // We might need to consider the  BroadCast case, where the same data is send to all the remote of the object
        // KM is the key_manager
        template < typename ITER, typename KM >
        void accumulateKeysOwnerScatter( ITER itb, ITER ite, KM & key_manager );
        // Question : should Owner send to himself ... not really send but everything done as if it had send ...
        //            This would insure same semantic as an MPI_BCast and perhaps simplify some algo

        /// each key in the iterator range send to all it's owners's remote key if the key is not owned
        /// Collective on the MPI_Comm of the exchanger
        /// Only avalaible for exchangeInformationPolicy with trait send_recv_keys_association_trait
        template < typename ITER, typename KM >
        void accumulateKeysOwnerGather( ITER itb, ITER ite, KM  & key_manager );
        // Question : should Owner send to himself ... not really send but everything done as if it had send ...
        //            This would insure same semantic as an MPI_BCast and perhaps simplify some algo

        /// each key in the iterator range send to all it's remote key
        /// Collective on the MPI_Comm of the exchanger
        template < typename ITER, typename KM >
        void accumulateKeysAllGather( ITER itb, ITER ite, KM & key_manager );
        /// each local object in the iterator range send to all it's remote object if the object is owned by the current proc. the data send to the remote object is not necesseraly the same.

        void clearKeys( );
        MPI_Comm getComm();

    private:
        typedef std::map < int, std::set < information_key_t > >                    info_key_recv_container_t;
        typedef std::map < int, std::map < information_key_t, information_key_t > > info_key_send_container_t;
        info_key_send_container_t send_keys_to;
        info_key_recv_container_t recv_keys_from;
        const MPI_Comm mpi_comm;
        int mpi_rank, mpi_size;

        template < typename KC, typename IM, typename COMMPATTERN, typename DATASTYLE >
        friend class exchangeInformationPolicy;

};

/// Container  that store the Keys in a send and receive keys association manner but only for a "from to" condition. Template on K, the type of key.
/// Self "from to" condition (i.e. from==to) is forbidden
template < typename K >
class xKeyContainerSendAndRecvFromTo
{
    public:
        typedef  K information_key_t;
        typedef send_and_recv_keys_association_trait keys_association_trait;

        xKeyContainerSendAndRecvFromTo( int from_, int to_, MPI_Comm _mpi_comm);

        /// each key in the iterator range that match the "from to" condition send/receive to/from "to"/"from" remote key
        /// Collective on the MPI_Comm of the exchanger
        template < typename ITER, typename KM >
        void accumulateKeys( ITER itb, ITER ite, KM & key_manager );

        void clearKeys( );
        MPI_Comm getComm();

    private:
        // Note: having map is unnecessary as send will always be for "to" and recv for "from" but
        // it stick to design of xKeyContainerSendAndRecv. This permit to use this class with exchangeInformationPolicy
        // having the benefit of already encoded communication pattern
        typedef std::map < int, std::set < information_key_t > >                    info_key_recv_container_t;
        typedef std::map < int, std::map < information_key_t, information_key_t > > info_key_send_container_t;
        info_key_send_container_t send_keys_to;
        info_key_recv_container_t recv_keys_from;
        const MPI_Comm mpi_comm;
        int from, to;
        int mpi_rank, mpi_size;

        template < typename KC, typename IM, typename COMMPATTERN, typename DATASTYLE >
        friend class exchangeInformationPolicy;

};

template < typename K,
           typename HASHFCN = std::hash < K >,
           typename EQUALKEY = std::equal_to < K >
           >
/// Container  that store the Keys in a send or recv keys association manner. Template on K, the type of key
//! note: 2 extras template parameter are required. Both got default. These are mandatory when K is not hashable by default hash
//! function. Same hold for EqualKey.
class xKeyContainerSendOrRecv
{
    public:
        typedef  K information_key_t;
        typedef send_or_recv_keys_association_trait keys_association_trait;


        xKeyContainerSendOrRecv( MPI_Comm _mpi_comm);

        /// Collective on the MPI_Comm of the exchanger
        template < typename ITER, typename KM >
        void accumulateKeys(ITER itb, ITER ite, KM & key_manager);
        template < typename OBJECTYPE, typename KM >
        void accumulateKeys(const OBJECTYPE & lo, KM & key_manager);

        void clearKeys( );
        MPI_Comm getComm();

    private:
        typedef std::map < int, std::unordered_set < information_key_t,HASHFCN,EQUALKEY > >                    info_key_container_t;
        info_key_container_t keys_from_to;
        const MPI_Comm mpi_comm;
        int mpi_rank, mpi_size;

        template < typename OBJECTYPE, typename KM >
        void accumulateKeysPrivate(const OBJECTYPE & lo, KM & key_manager);

        template < typename C, typename IM, typename COMMPATTERN, typename DATASTYLE >
        friend class exchangeInformationPolicy;
};

/// Class to pack information in a heterogeneous  buffer using MPI norm in underlying implementation
class xMpiInputBuffer
{
    public:
        inline xMpiInputBuffer(const MPI_Comm & _comm, int reserve = 0);
        inline int   pack(const void *inbuff, int incount, MPI_Datatype datatype);
        inline int   size() const;
        inline void *data();
    private:
        std::vector < char > dat;
        const MPI_Comm comm;
        mutable int position = 0;
};

/// Class to unpack information from a heterogeneous  buffer using MPI norm in underlying implementation
class xMpiOutputBuffer
{
    public:
        inline xMpiOutputBuffer(const MPI_Comm & _comm, int size);
        inline int         unPack(void *outbuff, int outcount, MPI_Datatype datatype) const;
        inline const void *data() const;
        inline void       *data();
        inline bool exhausted() const {return position >= (int) ( dat.size()); }
    private:
        const MPI_Comm comm;
        std::vector < char > dat;
        mutable int position = 0;
};

// exception
class xDataExchangerException : public std::exception
{
    public:
        xDataExchangerException(std::string message, std::string file, int line, std::string date, std::string time );
        ~xDataExchangerException() throw( ) override;
        const char * what() const throw( ) override;

    private:
        std::string msg;
};



} // end namespace

// Implementation
#include "xDataExchanger_imp.h"

#endif
