/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef XCOMM_H
#define XCOMM_H

#include <vector>
#include <cstddef>

namespace distmesh
{

/// Basic conector : give adress on a remote process memory space of a duplicate entity
//! This is a little container class that do nothing and is use as a bases of everything
//! It containe the adress and the id of the remote process
class xCommConector
{
    public:
        /// Default constructor wich set to zero adress and remote pid 
        xCommConector(void);
        /// Constructor which set  adress and remote pid  according to given parameter
        xCommConector(size_t rpid_,void *adress_);
        /// Methode to set  adress and remote pid  to values given as argument
        void set(size_t rpid,void *adress_);
        /// Methode to get  adress and remote pid  related to duplicate entity :
        //! \param [out] rpid id of the remote process
        //! \return Adress in remote process memory space of duplicated entity
        void *get(size_t & rpid);
    private:
        void *adress;
        size_t remote_pid;

};

/// Container of all duplicate entity of an entity on other proc
//! It stores some how instance of xCommConector that can be retrieved  by use of [] operator.
//! Underlyng storage is supposed to be continuous in memory an start at zero index up to 
//! size()-1 index. All xCommConector stored in this container are sorted in assending order of remote pid.
class xCommEntity
{
    public:
        ~xCommEntity();
        // public methodes ///////////////
        /// Methode to add a xCommConector to the collection of information stored in this instance
        void add(xCommConector conector);
        /// Methode to remove a xCommConector coresponding to remote proc id rproc from the collection of information stored in this instance
        //! Do nothing if no rproc present in the collection;
        //! This may transform this collection into an empty one.
        void remove(size_t rproc);
        /// Methode to retrive the actual size of the collection of xCommConnector stored in the instance
        size_t size(void);
        /// Methode to retrive the the \f$*idx^{th}\f$ xCommConnector stored in the instance
        xCommConector & operator [] (size_t idx);
        /// Methode to see if  first remote proc id is greater than proc_id argument 
        bool firstRemotIsGreaterThen(size_t proc_id);
        // public members ////////////////
    private:
        std::vector < xCommConector > connections;
};


} // end namspace
#endif
