/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef _XDISTVECTOR_H
#define _XDISTVECTOR_H

// std
#include <iostream>
#include <vector>
#include <cmath>

// xlinalg
#include "xDistIndex.h"

// distmesh
#include "xPartitionManager.h"
#include "xDataExchanger.h"

namespace xlinalg
{
/*!
 * This class correspond to a distributed vector concept. A vector (called V hereafter) of a certain size is dispatched on some proc.
 * All proc may have disconnect or overlapping part of this vector.The distributed nature of this vector is given at construction
 * time by a xDistIndex instance (called dist_index_ hereafter). This instance gives also the communicator on which this vector acte.
 * Part of the vector controlled by one process is compacted locally in an array corresponding to packed index order of dist_index_.
 *
 * Information stored in this vector are controlled by two states :
 *           - The global state correspond to a situation where a term present on a process reflect the value of this term in V.
 *           Global state is useful when you need to do matrix vector product with unassembled overlapping distributed matrix.
 *           With Xfem it is also needed to write or read, to or from, distributed double manager.
 *           - The local state correspond to a situation where local/owned (see xDistIndex) terms reflect the value of those term in V
 *           and all other term represent a contribution of those term in V. What is insured though is that the sum of all these
 *           overlapping contribution from all proc gives the term of V. This state is the natural state when inserting values in
 *           the vector. It is also mandatory for some algebraic operation like dot product where overlapping terms have to be
 *           "reduced" to be able to do correct local dot product.
 *
 *  In local state the constness have a special interpretation. Local term may be modified by some method but the sum remain
 *  the same.
 *
 *  Having those two state in mind one will understood that an other state must be added to control what can be done one the
 *  vector. Says inserting new value locally for example is clearly incompatible with global state.
 *  This extras state that will be call insert mode hereafter is controlled by user. If it is "on" user can modify/insert values
 *  in/to the vector but he can't do algebraic operation with this vector. By turn it "off" user can't insert value anymore but
 *  all algebraic method are then available. For sake of performance those control on insert mode are done with assert. Say in
 *  realese version there is no more control of insert mode status in methods of the class. This supose that user validate it's
 *  implementation (when in turn on or off the insert mode) in debug mode before passing to realese (fast) compilation version.
 *
 *  Method that work with insert mode "on" :
 *  - AddVal
 *  - getVal (const)
 *  - getVal
 *  - operator [](const)
 *  - operator []
 *  - at (const)
 *  - at
 *  - data (const)
 *  - data
 *  - switchInsertModeOff
 *  - switchInsertModeOffFromUserGlobalSetting
 *  - isInInsertModeOn
 *  - begin
 *  - end
 *  - getDistIndex
 *
 *  Method that work with insert mode "off" :
 *  - operator =
 *  - getVal (const)
 *  - operator [](const)
 *  - at (const)
 *  - data (const)
 *  - switchToGlobalValue
 *  - switchToLocalValue
 *  - switchInsertModeOn
 *  - isInInsertModeOn
 *  - dot
 *  - axpy
 *  - scal
 *  - nrm2
 *  - operator +,-,*,+=,-=
 *  - getDistIndex
 *
 *  Method that work with insert mode "on" or "off" with possible auto commuting status:
 *  - scatter
 *  - gather
 *
 *  NOTE: This storage may be treated with holes in indexing which make it a sparse storage as those unspecified index
 *  correspond in algebra computation to zero value.
 *  But gather and scatter method are in default in this case. Null term, i.e. non present index, are not properly taken into
 *  account for now. For gather and scater, you may also have problem  in interpretation of V vector if your indexing is not
 *  in ascending order across increasing process id. For example gathered information will correspond to dist_index order but
 *  may not be in increasing order of V indexation.
 */


template<typename VT = double>
class xDistVector
{
    private:
        typedef std::vector < VT > Vector;
    public:
        typedef typename Vector::iterator iterator;
        typedef typename Vector::const_iterator const_iterator;
        typedef typename Vector::reference reference;
        typedef typename Vector::reference const_reference;
        typedef typename Vector::value_type value_type;


        /*! Unique explicit constructor
         * This method give an instance in insert mode "on" and local state by default.
         * All stored values are set to 0.0.
         *   \param [in] dist_index_  is giving the indexing structure to adopt with this xDistVector instance. It provide also the communicator on which  vector will act.
         *   Not that you must use the same xDistIndex instance for all vector that you want to use together for algebraic operations
         */
        xDistVector( const xDistIndex & dist_index_ );

        /*! Default copy constructor is modified to :
         *  - set correctly internal menber which must be individualy copied
         *  \param [in] other an instance of an other xDistVector.
         */
        xDistVector(const xDistVector<VT> & other);

        /*! This method is creating a new intance of xDistVector based on same xDistIndex as the one who call it. 
         * This new instance is dynamicaly allocated and returned as a pointer. It is user responsability to dealocate this memory.
         * The given instance is in insert mode "off" and local state by default. All stored value are set to 0.0.
         * It give a way to use the cloned instance in algebra direcly.
         *
         * This method may be called with insert mode "on" or "off".
         */
        xDistVector<VT>* clone() const;

        /*! This method copy distributed vector given in argument into current vector
         * This method have to be called with insert mode "off"
         * \param [in] rhs other vector to copy (in any state)
         * \return A reference to this object
         */
        inline xDistVector<VT> & operator=(const xDistVector<VT> & rhs);

        /*! Set vector V entries to a specifique value alpha.
         * It is not a collective operation ! User have to maintain coherence. (see comment of scal)
         *
         * This method have to be called with insert mode "off".
         * \param [in] alpha the fixed value
         * \return A reference to this object
         */
        xDistVector<VT> & operator=(const VT & alpha);


        /*! This method add a value val at global index i in vector V.
         * This method is mandatory to make this class available for xfem Assemble functionality.
         *
         * This method have to be called with insert mode "on".
         * \param [in] i global index in Fortran indexing
         * \param [in] val local value to add to ith term of V
         */
        inline void AddVal( const int & i, const VT & val);

        /*! This method get a const reference to value at global index i in vector V.
         * If vector is in global state this value is the real value of ith term of V.
         * If vector is in local state this value is only a local contribution of this ith term of V.
         * - In Release version throw an exception if index is not in this process (loose checking).
         * - In Debug version assert stop execution if index is not in this process (strong checking).
         *
         * This method may be called with insert mode "on" or "off".
         * \param [in] i global index in Fortran indexing
         * \return A const reference to the value at global index i in vector V.
         */
        inline const VT & getVal(xDistIndex::idx_t i) const;

        /*! This method get a reference to value at global index i in vector V.
         * This instance is in local state and this value is only a local contribution of this ith term of V.
         * - In Release version throw an exception if index is not in this process (loose checking).
         * - In Debug version assert stop execution if index is not in this process (strong checking).
         *
         * This method have to be called with insert mode "on".
         * \param [in] i global index in Fortran indexing
         * \return A reference to the value at global index i in vector V.
         */
        inline VT & getVal(xDistIndex::idx_t i);

        /*! This method get a const reference to value at local index i in local storage.
         * Be careful on what you do with this operator as i is meaningless. You 'll have to use dist_index_ to
         * find out what does i represents.
         * If vector is in global state this value is the real value of a term of V.
         * If vector is in local state this value is only a local contribution of a term of V.
         * - In Release and Debug version don't throw an exception nor stop on assert if i is out of range.
         *
         * This method may be called with insert mode "on" or "off".
         * \param [in] i local index in C indexing of local packed storage
         * \return A const reference to the  ith value in local packed storage.
         */
        inline const VT & operator[](xDistIndex::idx_t i) const;

        /*! This method get a reference to value at local index i in local storage.
         * Be careful on what you do with this operator as i is meaningless. You 'll have to use dist_index_ to
         * find out what does i represents.
         * This instance is in local state and this value is only a local contribution of a term of V.
         * - In Release and Debug version don't throw an exception nor stop on assert if i is out of range.
         *
         * This method have to be called with insert mode "on".
         * \param [in] i local index in C indexing of local packed storage
         * \return A const reference to the  ith value in local packed storage.
         */
        inline VT & operator[]( xDistIndex::idx_t i);

        /*! See operator[] const
         * - In Release and Debug version throw an exception if i is out of range.
         */
        inline const VT & at(xDistIndex::idx_t i) const;

        /*! See operator[]
         * - In Release and Debug version throw an exception if i is out of range.
         */
        inline VT & at(xDistIndex::idx_t i);

        /*! This method give a const pointer to begin of local storage.
         * Be careful on what you do with this operator. Having access to underlying private memory space
         * may be dangerous. Adding friend class to access to this memory would have more sens as some how more secure.
         *
         * This method may be called with insert mode "on" or "off".
         * \return A const pointer to local packed storage.
         */
        inline const VT *  data() const;

        /*! This method give a pointer to begin of local storage.
         * Be careful on what you do with this operator. Having access to underlying private memory space
         * may be dangerous. Adding friend class to access to this memory would have more sens as some how more secure.
         *
         * This method have to be called with insert mode "on".
         * \return A pointer to local packed storage.
         */
        inline VT *  data();


        /*! This method give a iterator of the begin of the local storage.
         *
         * This method have to be called with insert mode "on".
         * \return A iterator to begin of local packed storage.
         */
        inline iterator begin();

        /*! This method give a iterator of pass the end of the local storage.
         *
         * This method have to be called with insert mode "on".
         * \return A iterator to pass theend  of local packed storage.
         */
        inline iterator end();

        /*! This method give a const iterator of the begin of the local storage.
         *
         * This method have to be called with insert mode "on".
         * \return A const iterator to begin of local packed storage.
         */
        inline const_iterator begin() const;

        /*! This method give a const iterator of pass the end of the local storage.
         *
         * This method have to be called with insert mode "on".
         * \return A const iterator to pass theend  of local packed storage.
         */
        inline const_iterator end() const;



        /*! Reset vector V entries to zero.
         * It is a collective operation.
         *
         * This method may be called with insert mode "on" or "off".
         */
        inline void Zerodata();


        /*! Printing in out the contend of this instance
         * It is a local operation. It is is just a snapshot of the current value status.
         * Very helpful for debugging it might be replace in the future by something following a well know format
         *
         * This method may be called with insert mode "on" or "off".
         */
        inline void Printdata( std::ostream  & out) const;


        //algebra
        /*! Compute the sum of this vector with another (rhs) and return the result in a new instance
         * It is a collective operation which may involve communication (mainly when vectors are not in the same state) .
         *
         * This method have to be called with insert mode "off"
         * \return A new vector sum of the operand of the +
         */
        inline xDistVector<VT> operator +(const xDistVector<VT> &rhs) const ;

        /*! Compute the difference of this vector with another (rhs) and return the result in a new instance
         * It is a collective operation which may involve communication (mainly when vectors are not in the same state) .
         *
         * This method have to be called with insert mode "off"
         * \return A new vector difference of the operand of the -
         */
        inline xDistVector<VT> operator -(const xDistVector <VT>&rhs) const;

        /*! Compute the product (dot product) of this vector with another (rhs) and return the result
         * It is a collective operation which involve communication.
         *
         * This method have to be called with insert mode "off"
         * \return the dot product of the operand of the *
         */
        inline VT operator *(const xDistVector<VT> &rhs) const;

        /*! Compute the product of this vector by a value (alpha) and return the result in a new instance
         * see discussion on scal
         *
         * This method have to be called with insert mode "off"
         * \return A new vector scaled by alpha
         */
        inline xDistVector<VT> operator *(const VT &alpha) const;

        /*! Add to this vector another (rhs) and return  a reference to itself
         * It is a collective operation which may involve communication (mainly when vectors are not in the same state) .
         *
         * This method have to be called with insert mode "off"
         * \return this instance after operation
         */
        inline xDistVector<VT> & operator +=(const xDistVector<VT> &rhs);

        /*! Subtract to this vector another (rhs) and return  a reference to itself
         * It is a collective operation which may involve communication (mainly when vectors are not in the same state) .
         *
         * This method have to be called with insert mode "off"
         * \return this instance after operation
         */
        inline xDistVector<VT> & operator -=(const xDistVector<VT> &rhs);

        /*! Compute the norm 2 of this vector
         * It is a collective operation which involve communication.
         *
         * This method have to be called with insert mode "off"
         * \return The norm
         */
        inline double nrm2( ) const; // warn involve communication

        /*! Compute the sum of this vector with the one given in argument scaled by "a" argument.
         * It is a collective operation which may involve communication (mainly when vectors are not in the same state) .
         *
         * This method have to be called with insert mode "off" and its vector argument must be in the same state.
         * \param [in] a the scaling value
         * \param [in] x a vector to compute with
         * \return nothing
         */
        void axpy(const VT &a, const xDistVector<VT> & x);

        /*! Multiply vector by a value
         * It is not a collective operation but it is user responsibility to maintain coherence.
         *
         * Says if one proc is scaling its part and the vector is in local mode it is not problematic.
         * If vector is in global mode, if you do a local scale its wrong as global values are not maintained across all proc.
         * For now if user want to used it it's is responsibility if he don't use this method collectively to switch his
         * vector in local mode.
         * If user do a collective call there is no issues.
         *
         * Note: This choice correspond to the aim of not forcing collective context with a barrier for example and not forcing
         * passing in local mode if not needed (collective with global mode is safe)
         *
         * This method have to be called with insert mode "off"
         * \param [in] alpha the scaling value
         */
        void scal(const VT & alpha);


        /*! Compute the dot product of two vectors
         * It is a collective operation which involve communication.
         *
         * This method have to be called with insert mode "off" and its argument must be in the same state.
         * \param [in] y a vector to compute with
         * \return The dot product "this.y"
         */
        VT dot( const xDistVector<VT> &y) const;

        /*! Compute the component by component product of two vectors x and y and 
         * store the result in calling instance.
         * It is a collective operation which may involve communication.
         *
         * This method have to be called with insert mode "off" and its argument must be in the same state.
         * \param [in] x a vector to compute with. 
         * \param [in] y a vector to compute with. 
         */
        void componentProduct( const xDistVector<VT>& x, const xDistVector<VT> &y);

        /*! Gather all values to form global vector V on  target_proc processor
         * Vector V what ever it contains will be overwritten and eventually resized on target_proc.
         * If user want to keep V in same memory definition he must have set V to appropriate dimension
         * to avoid any resize call.
         *
         * This method may be called with insert mode "on" or "off" (auto commuting status to insert mode "off" during operation).
         * \param [in,out] V vector where values will be gathered (see NOTE on general description of the class)
         * \param [in] target_proc processe on which we gone gather V
         */
        void gather( std::vector < VT > & V, int target_proc = 0) const;

        /*! Scatter all values from a global vector V on source_proc  processor to all local packed container.
         * Vector V must be of the correct size otherwise a exception is raised.
         *
         * This method may be called with insert mode "on" or "off" (auto commuting status to insert mode "on" during operation).
         * \param [in] V vector from which values will be scattered (see NOTE on general description of the class)
         * \param [in] source_proc processe which have proper definition of V
         */
        void scatter( std::vector < VT > & V, int source_proc = 0);

        /*! Put vector in global status : all values for a given index are the same across all proces.
         * This may involve communication.
         *
         * This method have to be called with insert mode "off".
         */
        void switchToGlobalValue();

        /*! Put vector in local status : sum of all value for a given index is correct but localy value is just a part of it
         *
         * This method have to be called with insert mode "off".
         */
        void switchToLocalValue();

        /*! This method switch to insert mode "on"
         *
         * This method have to be called with insert mode "off".
         */
        void switchInsertModeOn();

        /*! This method switch to insert mode "off"
         *
         * This method have to be called with insert mode "on".
         */
        void switchInsertModeOff();

        /*! This method switch to insert mode "off" but consider that 
         * user set value localy with global value. I.e. local values not owned by this proc are
         * set to the same value as the one set in process owning it.
         *
         * Typical usage: Using some custom mechanism user compute global value and use then to 
         *   set xDistVector instance  in insert mode "on" (via non const begin/end iterator for example).
         *   Then by using this method user may turn this instance in insert mode "off" and switch status
         *   to GlobalValue without any comunication.
         *
         * It is user responsability to insure that he provides really a set of value coresponding to Global state.
         * If not all remaing operation will be wrong ....
         *
         * This method have to be called with insert mode "on".
         */
        void switchInsertModeOffFromUserGlobalSetting();

        /*! This method check status of insert mode and return true if vector is in insert mode "on"
         *
         * This method may be called with insert mode "on" or "off".
         * \return True is  vector is in insert mode "on"
         */
        bool isInInsertModeOn() const ;

        /*! Extra friend function to fullfil requirement of BasicIterSolver interface at least
         * Note: all instance ave to be in insert mode "off"
         */
        //! result = v1 + c2*v2
        template<typename UT>
        friend void add(const xDistVector<UT> & v1, const UT & c2, const xDistVector<UT> & v2, xDistVector<UT> & result);
        //! v1 = v1 + c2*v2
        template<typename UT>
        friend void add(xDistVector<UT> & v1, UT c2, const xDistVector<UT> & v2 );
        //! result = c1*v1 + c2*v2
        template<typename UT>
        friend void add(const UT & c1, const xDistVector<UT> & v1, const UT & c2, const xDistVector<UT> & v2, xDistVector<UT> & result);
        //! result = alpha(v1 + v2)
        template<typename UT>
        friend void add(const UT & alpha, const xDistVector<UT> & v1, const xDistVector<UT> & v2, xDistVector<UT> & result);
        //! result = v1 + v2 + v3
        template<typename UT>
        friend void add(const xDistVector<UT> & v1, const xDistVector<UT> & v2, const xDistVector<UT> & v3, xDistVector<UT> & result);
        //! result = v1 - v2
        template<typename UT>
        friend void subtract(const xDistVector<UT> & v1, const xDistVector<UT> & v2, xDistVector<UT> & result);
        //! return the inner product of v1 and v2
        template<typename UT>
        friend UT InnerProduct(const xDistVector<UT> & v1, const xDistVector<UT> & v2);

        /*! Get const reference to the xDistIndex instance used to create
         * distributed vector
         *
         * This method may be called with insert mode "on" or "off".
         */
        const xDistIndex& getDistIndex() const;
    private:
        // private types
        class xKeyManagerDistVector
        {
            public:
                // mandatory traits
                typedef VT * information_key_t;

                 xKeyManagerDistVector(const xDistIndex &dist_index_,Vector &data_);

                // mandatory methods
                information_key_t localObjectKey( const xDistIndex::idx_t & lo);
                xtool::xConstPartitionObject < xDistIndex::idx_t > getConstPartitionObject( const xDistIndex::idx_t & lo);
                information_key_t remoteObjectKey(const xtool::xRemoteObject < xDistIndex::idx_t > & ro, const xDistIndex::idx_t & lo);
            private:
                const xDistIndex &dist_index;
                Vector &dat;


        };

        class xInfoManagerDistVectorReduce
        {
            public:
                typedef xtool::homogeneous_data_style_trait data_style_trait;
                typedef xtool::send_and_recv_keys_communication_trait communication_trait;
                typedef typename xKeyManagerDistVector::information_key_t information_key_t;
                typedef VT  information_t;
                information_t getInfo(information_key_t key, int sendto);
                void setInfo(information_key_t key, const information_t &info, int receivedfrom);
        };

        class xInfoManagerDistVectorDoubleReduce 
        {
            public:
                typedef xtool::homogeneous_data_style_trait data_style_trait;
                typedef xtool::send_and_recv_keys_communication_trait communication_trait;
                typedef typename xKeyManagerDistVector::information_key_t information_key_t;
                typedef std::array<VT, 2 >  information_t;
                xInfoManagerDistVectorDoubleReduce (Vector &data_, Vector &other);
                information_t getInfo(information_key_t key, int sendto);
                void setInfo(information_key_t key, const information_t &info, int receivedfrom);
            private:
                Vector &dat;
                Vector &other;
        };

        class xInfoManagerDistVectorSet
        {
            public:
                typedef xtool::homogeneous_data_style_trait data_style_trait;
                typedef xtool::send_and_recv_keys_communication_trait communication_trait;
                typedef typename xKeyManagerDistVector::information_key_t information_key_t;
                typedef VT  information_t;
                information_t getInfo(information_key_t key, int sendto);
                void setInfo(information_key_t key, const information_t &info, int receivedfrom);
        };
        enum status_t {RESET = 0,REDUCED = 1,GLOBAL = 2,FLAGRG = 3, INSERTOFF = 4,UNKNOWN};
        // private member dat
        Vector dat;
        const xDistIndex &dist_index;
        xInfoManagerDistVectorReduce reduce_info;
        xInfoManagerDistVectorSet set_info;
        xtool::xKeyContainerSendAndRecv < typename xKeyManagerDistVector::information_key_t >  reduce_keys;
        xtool::xKeyContainerSendAndRecv < typename xKeyManagerDistVector::information_key_t >     set_keys;
        char status;
        int proc_id;
        // private member function
        void reduceOnOwner();
        void doubleReduceOnOwner(xDistVector<VT> &other);
};

} // end of namespace

#include "xDistVector_imp.h"

#endif







