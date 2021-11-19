/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/
#ifndef _XDISTBLOCKDENSEMATRIX_H
#define _XDISTBLOCKDENSEMATRIX_H

// std
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// xlinalg
#include "xDistIndex.h"
#include "xDistVector.h"
#include "xTraitsMatrix.h"

// xTool
#include "xDataType.h"
#include "xPartitionManager.h"
#include "xDataExchanger.h"
#include "xMPIDataType.h"

namespace xlinalg
{
/*!
 * This class correspond to a distributed block dense square matrix concept. This may be a full dense matrix or have extra diagonal block that
 * are null. Each process have square dense block corresponding to part of the global matrix.
 * All proc may have disconnect or overlapping part of this matrix.The distributed nature of this matrix is given at construction
 * time by a xDistIndex instance (called dist_index_ hereafter). This instance gives also the communicator on which this matrix act.
 * Part of the matrix controlled by one process is compacted locally in an dense matrix corresponding to packed index order of dist_index_
 * in both dimension. This imply that this local matrix is square.
 *
 *  NOTE: This storage may be treated with holes in indexing which make it a sparse storage (block nature was also doing so)
 *  as those unspecified index correspond in algebra computation to zero value.
 *
 */
template < typename T,
           typename DEFINED,
           typename PATTERN
           >
class xDistBlockDenseMatrix
{
    private:
        typedef std::vector < T > Vector;
    public:
        // defined by the user
        typedef T matrix_data_type;
        typedef PATTERN matrix_pattern;
        typedef DEFINED matrix_defined; // maybe for future use with Lapack
        //
        // fixed
        typedef xTraitMatrixDense matrix_storage;
        typedef xTraitMatrixFindex matrix_indexing;
        typedef xTraitMatrixForcedAssemblyOnZero matrix_assembly_on_zero;
        typedef xDistIndex::idx_t idx_t;
        typedef typename Vector::iterator iterator;
        typedef typename Vector::const_iterator const_iterator;
        typedef typename Vector::reference reference;
        typedef typename Vector::reference const_reference;
        typedef typename Vector::value_type value_type;


        /*! Unique explicit constructor
         * This method give an instance where all stored values are set to 0.0.
         *   \param [in] dist_index_  is giving the indexing structure to adopt with this xDistBlockDenseMatrix instance.
         *                            It provide also the communicator on which  vector will act.
         *
         *   Note that you must use the same xDistIndex instance for all vector that you want to use with this matrix for algebraic
         *   operations
         *   Note that your local dense block stored in this instance is supposed to be square of dimension n_bloc x n_bloc where
         *   n_bloc is the size of the packed index defined in dist_index_. The size of the global matrix is n_glob x n_glob where
         *   n_glob is the size of all index of dist_index_.
         *   For symmetric matrix the size of the local dense is the same but only lower triangular part is used
         */
        xDistBlockDenseMatrix( const xDistIndex & dist_index_ );

        /*! Method to set  real values of matrix term being owned by the current process. By real we mean here the sum of all contribution 
         * from remote counterpart( in sens of xDistIndex). Remote counterpart terms are set to zero after this call so that sum across all
         * process of those terms remain the same as before the call.
         */
        void reduceOnOwner();

        /*! Printing in out the contend of this instance
         * It is a local operation. It is is just a snapshot of the current value status.
         * Very helpful for debugging it might be replace in the future by something following a well know format
         *
         */
        inline void Printdata( std::ostream  & out) const;

        /*! method needed to connect this matrix to a solver or whatever
         *
         *   \param [out] ng  dimension of the global matrix
         *   \param [out] nl  dimension of the local dense block
         *   \param [out] data  pointer to the beginning of the storage of the dense block in column major ordering (i+nl*j give index of i,j term)
         */
        void getMemoryAccess(int &ng,int &nl,T **data);

        /// do Y = alpha.A.X + beta.Y  where A is this matrix
        //! For now this use dense algorithm strategy without using current status (if reduced we are multiplying by potentially many zero terms)
        //
        void gemv( T alpha, T beta,const xlinalg::xDistVector<T> & X, xlinalg::xDistVector<T> & Y);


    private:
        // private types
        typedef std::map < const idx_t *, T * >  index_column_asso_t;
        class pairIntHash
        {
            public:
                std::size_t operator()(const std::pair < int, int > &x) const
                {
                    return string_hash( std::to_string(x.first)+" "+std::to_string(x.second));
                }
            private:
                std::hash < std::string > string_hash;
        };
        typedef std::unordered_map < std::pair < int,int >,index_column_asso_t,pairIntHash > index_column_asso_container_t;
        class xKeyManagerDistBlockDenseMatrix
        {
            public:
                // data type: type of object used to discribe keys
                typedef idx_t data_t;

                // mandatory traits
                typedef T * information_key_t;

                xKeyManagerDistBlockDenseMatrix(const xDistIndex &dist_index_,Vector &data_);

                // mandatory methods
                information_key_t localObjectKey( const data_t & lo);
                xtool::xConstPartitionObject < xDistIndex::idx_t > getConstPartitionObject( const data_t & lo);
                information_key_t remoteObjectKey(const xtool::xRemoteObject < xDistIndex::idx_t > & ro, const data_t & lo);
            private:
                const xDistIndex &dist_index;
                Vector &dat;

        };
        class xInfoManagerDistBlockDenseMatrixReduce
        {
            public:
                xInfoManagerDistBlockDenseMatrixReduce(Vector &dat_,index_column_asso_container_t &exchanged_column_terms_);
                typedef xtool::nonhomogeneous_data_style_trait data_style_trait;
                typedef xtool::send_and_recv_keys_communication_trait communication_trait;
                typedef typename xKeyManagerDistBlockDenseMatrix::information_key_t information_key_t;
                void getInfo(information_key_t key, xtool::xMpiInputBuffer & buff,int sendto);
                void setInfo(information_key_t key,  const xtool::xMpiOutputBuffer & buff, int receivedfrom);
                size_t getApproxDataSize(void);
                void setMoy(int moy_);
            private:
                Vector &dat;
                index_column_asso_container_t &  exchanged_column_terms;
                T zero;
                MPI_Datatype mpi_data_t;
                int moy;

        };
        class xInfoManagerDistBlockDenseMatrixSet
        {
            public:
                xInfoManagerDistBlockDenseMatrixSet(Vector &dat_,index_column_asso_container_t &exchanged_column_terms_);
                typedef xtool::homogeneous_data_style_trait data_style_trait;
                typedef xtool::send_and_recv_keys_communication_trait communication_trait;
                typedef typename xKeyManagerDistBlockDenseMatrix::information_key_t information_key_t;
                void getInfo(information_key_t key, xtool::xMpiInputBuffer & buff,int sendto);
                void setInfo(information_key_t key,  const xtool::xMpiOutputBuffer & buff, int receivedfrom);
                size_t getApproxDataSize(void);
                void setMoy(int moy_);
            private:
                Vector &dat;
                index_column_asso_container_t &  exchanged_column_terms;
                MPI_Datatype mpi_data_t;
                int moy;
        };

        enum status_t {RESET = 0,REDUCED = 1,GLOBAL = 2,FLAGRG = 3, INSERTOFF = 4,UNKNOWN};
        // private member
        int n_bloc;
        Vector dat;
        Vector state_boundary_dat;
        index_column_asso_container_t exchanged_column_terms;
        const xDistIndex &dist_index;
        xInfoManagerDistBlockDenseMatrixReduce reduce_info;
        xInfoManagerDistBlockDenseMatrixSet set_info;
        xtool::xKeyContainerSendAndRecv < typename xKeyManagerDistBlockDenseMatrix::information_key_t >  reduce_keys;
        xtool::xKeyContainerSendAndRecv < typename xKeyManagerDistBlockDenseMatrix::information_key_t >     set_keys;
        char status;
        int proc_id;
        // private member function
};

//===== xPolicyxDistBlockDenseMatrixGemv ================================================================================
template < typename T, typename PATTERN >
class xPolicyxDistBlockDenseMatrixGemv
{
    public:
        static void  mv( const int & n, const double & alpha, const double*A, const int & LDA, const double *X, const double & beta, double *Y);
};

} // end of namespace

#include "xDistBlockDenseMatrix_imp.h"

#endif







