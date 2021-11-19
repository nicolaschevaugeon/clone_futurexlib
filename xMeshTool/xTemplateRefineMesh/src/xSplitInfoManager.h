/* 
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE & LICENSE files for terms 
   and conditions.
*/

#ifndef XSPLITINFOMANAGER_H
#define XSPLITINFOMANAGER_H
#include <set>
#include <array>



//xfem::DistMesh
#include "xDataExchangerTraits.h"


namespace xmeshtool
{

  template < typename PM, typename T >
    class splitKeyManagerSendAndReceive
  {
  public:
    /// basic constructor which retrive partition manager associated to this class
    splitKeyManagerSendAndReceive( PM &partman_ );

    // mandatory types for information class familly
    typedef const T * information_key_t;

    // mandary method for key class familly
    information_key_t localObjectKey( const T & o);
    information_key_t remoteObjectKey(const xtool::xRemoteObject < T > & ro, const T  &lo );
    xtool::xConstPartitionObject < T > getConstPartitionObject(const T &e);
  protected:
    PM &partman;
  };
  template < typename PM, typename T >
    class splitKeyManagerSendOnly
  {
  public:
    /// basic constructor which retrive partition manager associated to this class
    splitKeyManagerSendOnly( PM &partman_ );

    // mandatory types for information class familly
    typedef const T * information_key_t;


    // mandary method for key class familly
    information_key_t localObjectKey( const T & o);
    std::set < int > getMessageRanks( const information_key_t & lk);
  protected:
    PM &partman;
  };

  /// One Level information manager
  //! Exchange on proc boundary the fact that, for a given frontier, related elements in remot proc have to be split
  //! Send adress of remote object so that it is added to to_be_added on remote proc
  //! T is the type of the elements
  template < typename PM, typename T >
    class splitInfoManagerOneLevelInformation
  {
  public:
    /// Constructor
    splitInfoManagerOneLevelInformation(const PM &part_man_,std::set < T * > &to_be_added_);

    // mandatory types for information class familly
    typedef xtool::homogeneous_data_style_trait data_style_trait;
    typedef xtool::send_only_keys_communication_trait communication_trait;
    typedef const T * information_key_t;
    typedef T * information_t;

    // mandatory method for information class familly
    information_t getInfo(information_key_t key, int sendto);
    void setInfo(const std::vector < information_t > &info, int receivedfrom );
  protected:
    const PM &partman;
    std::set < T * > & to_be_added;
  };

  /// Treat new hanging edge/face information with treat_new_hanging_
  template < typename PM, typename T >
    class splitInfoManagerNewHangingTreatement
  {
  public:
    /// Constructor
    splitInfoManagerNewHangingTreatement(const PM &part_man_,std::function < void(T &) > &treat_new_hanging_);

    // mandatory types for information class familly
    typedef xtool::homogeneous_data_style_trait data_style_trait;
    typedef xtool::send_only_keys_communication_trait communication_trait;
    typedef const T * information_key_t;
    typedef T * information_t;

    // mandatory method for information class familly
    information_t getInfo(information_key_t key, int sendto);
    void setInfo(const std::vector < information_t > &info, int receivedfrom );
  protected:
    const PM &partman;
    std::function < void (T &) > treat_new_hanging;
  };

  /// update partition manager to take into account new hanging information (edges/faces) related to object (edge/faces) selected during accumulation
  //! It treat also split face for terminal element
  template < typename PM, typename T, typename DM, typename DM1, typename DM2 >
    class splitInfoManagerUpdatePMFE
  {
  public:
    /// Constructor
    splitInfoManagerUpdatePMFE( PM &part_man_,DM & hd, DM & hu, DM1 &st1, DM2 &st2);

    // mandatory types for information class familly
    typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
    typedef xtool::send_and_recv_keys_communication_trait communication_trait;
    typedef const T * information_key_t;


    // mandatory method for information class familly
    void getInfo(information_key_t key, xtool::xMpiInputBuffer & buff, int sendto);
    void setInfo(information_key_t key, const xtool::xMpiOutputBuffer & buff, int receivedfrom);
    size_t getApproxDataSize();

  protected:
    PM & partman;
    DM & hanging_down;
    DM & hanging_up;
    DM1 & sub_t1;
    DM2 & sub_t2;
  };

  /// Count conection to identify dangling entity
  template < typename T, template < typename >  class DM >
    class splitInfoManagerDanglingCount 
    {
    public:
      // mandatory types for information class familly
      typedef xtool::homogeneous_data_style_trait data_style_trait;
      typedef xtool::send_and_recv_keys_communication_trait communication_trait;
      typedef const T * information_key_t;
      typedef unsigned short int information_t;

      // specific to this class
      void setNbConnect(const information_key_t & k, const information_t &inf);
      void getNbConnect(const information_key_t & k, information_t &inf);
      void clearNbConnect(const information_key_t & k);


      // mandatory method for information class familly
      information_t getInfo(information_key_t key, int sendto);
      void setInfo(information_key_t key, const information_t &info, int receivedfrom);
    protected:
      DM < information_t > nb_component;
    };


} // end namspace

#include "xSplitInfoManager_imp.h"

#endif
