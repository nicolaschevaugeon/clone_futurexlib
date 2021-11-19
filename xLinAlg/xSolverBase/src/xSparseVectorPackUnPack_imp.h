/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifdef XSPARSEVECTORPACKUNPACK_H
#ifndef XSPARSEVECTORPACKUNPACK_IMP_H
#define XSPARSEVECTORPACKUNPACK_IMP_H

template <typename T>
void pack(const xlinalg::xSparseVector<T> &vect, xtool::xMpiInputBuffer &buff)
{
    int nb=vect.nnz();
    buff.pack(&nb,1,MPI_INT);
    if (nb)
    {
        buff.pack(&(*vect.beginIdx()),nb,MPI_INT);
        buff.pack(&(*vect.beginVal()),nb,xtool::xMPIDataType<T>());
    }
}
template <typename T>
void unPack(xlinalg::xSparseVector<T> &vect, const xtool::xMpiOutputBuffer &buff)
{
    int nb;
    buff.unPack(&nb,1,MPI_INT);
    if (nb)
    {
        std::vector<int> idx(nb);
        buff.unPack(idx.data(),nb,MPI_INT);
        std::vector<T> val(nb);
        buff.unPack(val.data(),nb,xtool::xMPIDataType<T>());
        vect.resetBy(idx.begin(),idx.end(),val.begin(),val.end());
    }
    else
        vect.clear();

}

#endif
#endif
