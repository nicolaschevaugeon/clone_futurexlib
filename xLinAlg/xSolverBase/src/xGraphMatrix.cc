/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <algorithm>
#include <cassert>
#include "xGraphMatrix.h"
#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>

using namespace std;

namespace xlinalg
{

xGraphMatrix::xGraphMatrix(int n_, int m_, bool sym_) : n(n_), m(m_), nnz(0), sym(sym_)
{
    // set vector size from number of colonnes
    matrix_struct.resize(m);

    // check :
    //
    // symetrique pattern used to fill graph (call AssembleGraph with one field)
    // -------------------------------------------------------------------------
    //
    // if sym=1 we store only half of the symetric pattern involved here
    // but if n#m we are facing 2 situations :
    //  n>m : it have sens, has we store matrix structure in sparse colone format,
    //        lower half representation is correct for the bloc 1:m x 1:m. The actual "field"
    //        is included in this sub bloc.
    //        Terms of block m+1:n x 1:m may exist from a other call to AssembleGraph and they remain untouched
    //  n<m : It doesn't have much sens has in this case terms of block 1:n x n+1:m vanish by
    //        the simple fact that lower bloc n+1:m x 1:n doesn't exist
    //
    //
    // unsymetrique pattern used to fill graph (call AssembleGraph with two different field)
    // -------------------------------------------------------------------------------------
    //
    // if sym=1 we are in the case where this unsymetric pattern involved here is stored in the
    //   lower part of a global symetric matrix
    //  if n#m we are facing 2 situations :
    //  n>m : it have sens, has we store matrix structure in sparse colone format,
    //        lower half representation is correct for the bloc 1:m x 1:m. The actual unsymetric
    //        pattern may have somme term in this sub bloc.
    //        Terms of block m+1:n x 1:m exist an the actual unsymetric pattern
    //        may have somme term in this sub bloc.
    //  n<m : It doesn't have much sens has in this case terms of block 1:n x n+1:m vanish by
    //        the simple fact that lower bloc n+1:m x 1:n doesn't exist
    //

    //
    if (sym && n < m)
    {
        std::cout << "In file "<< __FILE__ << " line " << __LINE__ << " compiled "<<__DATE__<<" at "<<__TIME__<<std::endl;
        std::cout << "Asking for symetric storage with n("<<n<<") < m("<<m<<") is not possible"<<std::endl;
        throw;
    }

    // setting addLinesFunction
    if (sym)
        addLinesFunction = std::bind(&xGraphMatrix::addLinesSym,this,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4);
    else
        addLinesFunction = std::bind(&xGraphMatrix::addLinesUnSym,this,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4);


}

// ReoderGraphFromReordering implementation --------------------------------------------------------------------------------------------
void xGraphMatrix::ReoderGraphFromReordering(int nr, int *perm_r, int nc, int *perm_c )
{
    // check
    if (nr != n || nc != m)
    {
        std::cout << "In file "<< __FILE__ << " line " << __LINE__ << " compiled "<<__DATE__<<" at "<<__TIME__<<std::endl;
        std::cout << "Dimensional problem ! You give for permutation vector the folowing sizes :\nnr="<<nr;
        std::cout << " nc="<<nc<<"\n but dimension of the matrix given so far is :\nn="<<n;
        std::cout << " m="<<m<<std::endl;
        throw;
    }
    // local
    int k,l;

    // Reorder ptr and reorder indices first pass : apply permutation on index
    vector < bool > flag(n,1);
    for (k = 0; k < n; ++k)
    {
        if (flag[k])
        {
            std::vector < vector < int > >::iterator itk = matrix_struct.begin();
            itk += k;
            l = perm_c[k];
            while (l != k)
            {
                std::vector < int >::iterator it0 = itk->begin();
                std::vector < int >::iterator it = it0;
                std::vector < int >::iterator itend = itk->end();
                for (; it != itend; ++it)
                {
                    ( *it ) = perm_r[( *it )];
                }
                sort(it0,itend);
                assert(flag[l]);
                std::vector < vector < int > >::iterator itl = matrix_struct.begin();
                itl += l;
                std::iter_swap(itk,itl);
                flag[l] = 0;
                l = perm_c[l];
            }
            std::vector < int >::iterator it = itk->begin();
            std::vector < int >::iterator itend = itk->end();
            for (; it != itend; ++it)
            {
                ( *it ) = perm_r[( *it )];
            }
            flag[k] = 0;
        }
    }

    // If sym=true a second pass to reoder indices is needed : re dispatch index in proper colonne as this
    // storage represent only half part of the matrix
    if (sym)
    {
        for (k = 0; k < n; ++k)
        {
            vector < int > &new_idx = matrix_struct[k];
            vector < int > old_idx;
            old_idx.swap(new_idx);
            std::vector < int >::iterator it = old_idx.begin();
            std::vector < int >::iterator itend = old_idx.end();
            for (; it != itend; ++it)
            {
                addSym(( *it ),k);
            }
        }
    }

    return;
}

// Add methode to the graphe  --------------------------------------------------------------------------------------------
void xGraphMatrix::addLinesSymBlock(const int il,const int nb_tot,const int j, int *line_idx0)
{
    addLinesFunction(il,nb_tot,j,line_idx0);
}
void xGraphMatrix::add(const int i,const int j)
{

    vector < int > &idx = matrix_struct[j];

    std::vector < int >::iterator itlower = lower_bound(idx.begin(),idx.end(),i);
    if (itlower == idx.end())
        idx.push_back(i);
    else if (( *itlower ) > i)
        idx.insert(itlower,i);


}
void xGraphMatrix::addSym(const int i,const int j)
{
    int k,l;
    if (i < j)
    {
        l = i;
        k = j;
    }
    else
    {
        k = i;
        l = j;
    }
    add(k,l);
}
void xGraphMatrix::addLinesUnsymBlock(const int nb,const int j, int *line_idx)
{
    int min = line_idx[0];

    vector < int > &idx = matrix_struct[j];
    const int ni = idx.size();

    if (ni)
    {
        //  test maximum
        std::vector < int >::iterator itend = idx.end();
        int max_idx = ( *( itend-1 ) );

        if (max_idx < min)
            idx.insert(itend,line_idx,line_idx+nb);
        else if (max_idx > min)
        {
            std::vector < int >::iterator it = idx.begin();
            vector < int > tmp;
            tmp.reserve(nb+ni);
            std::vector < int >::iterator itt = tmp.begin();
            itend = std::set_union(it,itend,line_idx,line_idx+nb,itt);
            idx.assign(itt,itend);
        }
        else
            idx.insert(itend,line_idx+1,line_idx+nb);
    }
    else
    {
        idx.resize(nb);
        memcpy ( (void *) ( &idx[0] ), (void *) line_idx, sizeof( int )*nb );
    }
    return;
}
void xGraphMatrix::addLinesUnSym(const int il,const int nb_tot,const int j, int *line_idx0)
{
    addLinesUnsymBlock(nb_tot,j,line_idx0);
}
void xGraphMatrix::addLinesSym(const int il,const int nb_tot,const int j, int *line_idx0)
{
    const int nb = nb_tot-il;
    int *line_idx = line_idx0+il;

    if (nb > 1)
    {
        addLinesUnsymBlock(nb,j,line_idx);
    }
    else
    {
        add(line_idx[0],j);
    }

    return;
}
// helper  --------------------------------------------------------------------------------------------
void xGraphMatrix::countNNZ()
{
    // reset to zero
    nnz = 0;

    // loop on column
    for (int i = 0; i < m; i++)
    {
        nnz += matrix_struct[i].size();
    }

}
void xGraphMatrix::remove(const int i,const int j)
{
    vector < int > &idx = matrix_struct[j];
    vector < int >::iterator it = find(idx.begin(), idx.end(), i);
    if (it != idx.end()) idx.erase(it);
}
int xGraphMatrix::getBandWidth()
{
    // local
    int d,l;

    if (sym)
    {
        int band = 0;
        // loop on column
        for (int j = 0; j < m; j++)
        {
            vector < int > &idx = matrix_struct[j];
            l = idx.size();
            if (l)
            {
                d = idx[l-1]-idx[0];
                if (d > band)
                    band = d;
            }
        }
        return ( 2*band+1 );
    }
    else
    {
        int band1 = 0;
        int band2 = 0;
        // loop on column
        for (int j = 0; j < m; j++)
        {
            vector < int > &idx = matrix_struct[j];
            l = idx.size();
            if (l)
            {
                std::vector < int >::iterator itlower = lower_bound(idx.begin(),idx.end(),j);
                if (itlower == idx.end())
                {
                    d = idx[l-1]-idx[0];
                    if (d > band1)
                        band1 = d;
                }
                else if (( *itlower ) > j)
                {
                    d = j-idx[0];
                    if (d > band1)
                        band1 = d;
                    d = idx[l-1]-( *itlower );
                    if (d > band2)
                        band2 = d;
                }
                else
                {
                    d = j-idx[0];
                    if (d > band1)
                        band1 = d;
                    d = idx[l-1]-j;
                    if (d > band2)
                        band2 = d;
                }
            }
        }
        return( band1+band2+1 );
    }
}
void xGraphMatrix::clear()
{

    // nota : A priori no memory leeks here as vector contained in matrix_struct
    //        are also clear in a more drastic way : destructor call
    matrix_struct.clear();

    n = m = nnz = 0;
}

void xGraphMatrix::print()
{
    int i,j;
    cout<<"Graph of a "<<n<<"x"<<m<<" matrix"<<endl;
    if (sym)
        cout<<"Graph is considered some how symetric and only lower part of it is stored"<<endl;
    else
        cout<<"Graph is unsymetrique and all terms are stored"<<endl;
    for (j = 0; j < m; ++j)
    {
        cout<<"Colonne "<<j<<" : "<<endl;
        vector < int > &idx = matrix_struct[j];
        const int ni = idx.size();
        if (ni)
        {

            cout<<ni<<" non zero termes :\n";
            for (i = 0; i < ni; ++i) cout<<"idx["<<i<<"]="<<idx[i]<<endl;
        }
        else
            cout<<"is empty"<<endl;
    }

}

void xGraphMatrix::printCOO(std::ofstream & os)
{
    int i,j;
    for (j = 0; j < m; ++j)
    {
        vector < int > &idx = matrix_struct[j];
        const int ni = idx.size();
        const int j1 = j+1;
        if (ni)
        {

            for (i = 0; i < ni; ++i) os<<(int) ( idx[i]+1 )<<" "<<j1<<" 1.\n";
        }
    }

}

// Permutation ----------------------------------------------------------------------------------------
void xGraphMatrix::getPermReverseCutHillMcKee(int nm, int *perm,int *invperm)
{
    if (n != m || !sym || nm != n )
    {
        std::cout << "In file "<< __FILE__ << " line " << __LINE__ << " compiled "<<__DATE__<<" at "<<__TIME__<<std::endl;
        if (n != m || !sym) std::cout<<"Only symetric square matrix is possible for getPermReverseCutHillMcKee. perm,invperm unchanged"<<endl;
        if ( nm != n )
        {
            std::cout << "Dimensional problem ! You give for permutation vector the folowing sizes :\nn="<<nm;
            std::cout <<"\n but dimension of the matrix given so far is :\nn="<<n;
            std::cout << " m="<<m<<std::endl;
            throw;
        }
        return;
    }

    // use boost implementation => transforme matrix_struct in boost graph
    //
    // nedeed types
    typedef boost::adjacency_list <
        boost::vecS, boost::vecS,boost::undirectedS,
        boost::property < boost::vertex_color_t, boost::default_color_type,
                          boost::property < boost::vertex_degree_t, int > > > Graph;
    typedef boost::graph_traits < Graph >::vertices_size_type vertices_index_t;
    //    typedef boost::graph_traits< Graph >::vertex_iterator      vertex_iterator_t;
    //typedef boost::graph_traits< Graph >::out_edge_iterator    edge_iterator_t;
    //typedef boost::graph_traits< Graph >::adjacency_iterator   adjacency_iterator_t;
    typedef boost::graph_traits < Graph >::vertex_descriptor gVertex;
    typedef boost::property_map < Graph, boost::vertex_index_t >::type property_map_type;
    //typedef boost::property_traits<boost::property_map< Graph, boost::vertex_index_t >::const_type>::value_type property_value_type;


    // local
    property_map_type old_index;
    vertices_index_t i,j,ni,nmt = (vertices_index_t) ( nm );



    // create reverse index permutation storage
    vector < gVertex > inv_perm(nmt);

    // contex to remove as soon as possible boost graph
    {
        // local
        Graph g(nmt);

        //fill graph
        for (j = 0; j < nmt; ++j)
        {
            vector < int > &idx = matrix_struct[j];
            if ((vertices_index_t) ( idx[0] ) != j)
            {
                std::cout << "In file "<< __FILE__ << " line " << __LINE__ << " compiled "<<__DATE__<<" at "<<__TIME__<<std::endl;
                std::cout << "First term of the matrix column should be the diagonal term !\n here it's ("<<idx[0]<<","<<j<<")\n";
            }
            ni = idx.size();
            for (i = 1; i < ni; ++i)
            {
                boost::add_edge(idx[i],j,g);
            }
        }

        // initial index of the vertice in the graph (property here are index (color))
        old_index = boost::get( boost::vertex_index, g );

        // reorder with cuthill mckee algo
        boost::cuthill_mckee_ordering( g, inv_perm.rbegin(),
                                       boost::get( boost::vertex_color, g ),
                                       boost::make_degree_map( g ) );

        // end context : graph destructed
    }

    // create permutation
    for (i = 0; i < nmt; ++i) invperm[ old_index[ inv_perm[i] ] ] = i;

    for (i = 0; i < nmt; ++i) perm[ invperm[ i ] ] = i;


}

} // end of namespace
