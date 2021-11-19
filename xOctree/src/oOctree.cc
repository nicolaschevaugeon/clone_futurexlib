/*
    octree is a subproject of  xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#include <cstdio>
#include <vector>
#include <cmath>
#include <iostream>
#include <iterator>
#include <algorithm>
#include "oOctree.h"
#include "oKeyManager.h"

using namespace std;

namespace xoctree {

void oOctree::construct ( const int *motif )
      {
    const bool debug = false;

    int powa=1;
    if(motif) powa = topo.powint(2, motif[0]+motif[1]+motif[2]);

    nsize[0] = powa;
    for (int l = 1; l <= LevelMax; ++l)
    {
      powa *= NbChildren;
      nsize[l] = powa;
      if (debug) cout << " size of level " << l << " is " << nsize[l] << endl;
    }
    
    ibegin[0] = 0;
    iend[0]   = nsize[0];
    for (int l = 1; l <= LevelMax; ++l)
    {
      ibegin[l] = iend[l-1];
      iend[l] = ibegin[l] + nsize[l];
      if (debug) cout << " for level " << l << " ibegin is " << ibegin[l] << " iend is " << iend[l] << endl;
    }
    
    CumulativeSize = iend[LevelMax];
    if (1) cout << "CumulativeSize " <<  CumulativeSize << endl;
    cells = new cell_type [ CumulativeSize ];
//     cells.resize(CumulativeSize);
//     cells = new std::vector<cell_type>[ CumulativeSize ];

    int imax = topo.index_max_for_levels[0][LevelMax];
    //cout<<"imax="<<imax<<endl;
    if(motif){
    // On impose qd meme le meme niveau max selon x, y et z
    int lmax0=*max_element(motif, motif+3);
    imax = topo.idx_max_for_level(LevelMax+lmax0,0) + lmax0;
    }

    
    ijk2_almost_octree = new int [imax + 1];
    
    int coefs_cart[topo.MAX_DEPTH];
    int size_cart, c;
    for (int i=0; i <= imax ; ++i)
    {
      topo.decompose (i, 2, coefs_cart, size_cart);
      topo.compose(c, NbChildren,  coefs_cart, size_cart);
      ijk2_almost_octree[i] = c;
    }

    //most refined level is active
    fill(begin(0),     end(LevelMax-1), 0);
    fill(begin(LevelMax), end(LevelMax), 1);
//throw;
  }

  oOctree::oOctree ( const oMapping& m_, int _LevelMax, int* per, bool construct_ )
    : mapping(m_), LevelMax (_LevelMax), Dim(m_.getDim()),
      NbChildren (topo.pow_base2[m_.getDim()]),  cellsSaved(0), topo(mapping.getTopo())
  {

    if(Dim==3) assert(LevelMax <= 10);
    if(Dim==2) assert(LevelMax <= 15);

    nsize= new int[topo.MAX_DEPTH];
    ibegin= new int[topo.MAX_DEPTH];
    iend= new int[topo.MAX_DEPTH];

    if (!per) fill(Periodicity, Periodicity+3, 0);
    else      copy(per, per+3, Periodicity);
    
    if(construct_) construct( );

  }

  oOctree::~oOctree ()
  {
    
    delete [] cells;
    delete [] ijk2_almost_octree;
    delete [] nsize;
    delete [] ibegin;
    delete [] iend;
  }


  void oOctree::activate  (oOctree::cell_type* cell, int level, const int* ijk)
  {
    iterator c = cell;
    int size = 1;
    *cell += 1;
    for (int l = level; l < LevelMax; ++l)
      {
        size *= NbChildren;
        c = getChildren(c, l);
        transform(c, c + size, c, bind2nd(std::plus<int>(), 1));
      }
    notify_derefinement(cell, level, ijk);
  }

  void oOctree::deactivate  (oOctree::cell_type* cell, int level, const int* ijk)
  {
    iterator c = cell;
    int size = 1;
    *cell -= 1;
    for (int l = level; l < LevelMax; ++l)
      {
        size *= NbChildren;
        c = getChildren(c, l);
        transform(c, c + size, c, bind2nd(std::minus<int>(), 1));
      }
    notify_refinement(cell, level, ijk);
  }



oOctree::cell_type* oOctree::getAncestor(int up, const oOctree::cell_type* cell, int level)
  {
    const bool debug = false;
    assert(level > 0);
    assert(up > 0);
    int Ioctree = cell - begin(level);
    oOctree::cell_type * ret = begin(level-up) + Ioctree/topo.generation_size[Dim-1][up];
    if (debug) cout << " the ancestor of Ioctree " << Ioctree
            << " at level " << level << " for " << up << " generation above is" <<
      Ioctree/topo.generation_size[Dim-1][up] << endl;

    return ret;
  }

  void oOctree::printActivity() const
  {
    for (int l = 0; l <= LevelMax; ++l)
      {
        cout << "activities at level " << l << " are " << endl;
        copy(begin(l), end(l), ostream_iterator<int>(cout, " "));
        cout << endl;
      }
  }


  void oOctree::viewState()
  {
    for (int l = 0; l <= LevelMax; ++l)
      {
        //       cout << "activities at level " << l << " are " << endl;
        copy(begin(l), end(l), ostream_iterator<int>(cout, " "));
        //       cout << endl;
      }
    cout << endl;
  }


  bool oOctree::is_on_boundary(int b, const int* ijk, int level) const
  {
    bool ret;
    if (Dim == 2)
      {
        if (b == 0)    ret =  (ijk[1] == 0);
      else if (b==1) ret =  (ijk[0] == topo.idx_max_for_level(level,0));
      else if (b==2) ret =  (ijk[1] == topo.idx_max_for_level(level,1));
      else if (b==3) ret =  (ijk[0] == 0);
      else throw -10;
    }
  else if (Dim == 3)
    {
      if (b == 0)    ret =  (ijk[1] == 0);
      else if (b==1) ret =  (ijk[0] == topo.idx_max_for_level(level,0));
      else if (b==2) ret =  (ijk[1] == topo.idx_max_for_level(level,1));
      else if (b==3) ret =  (ijk[0] == 0);
      else if (b==4) ret =  (ijk[2] == 0);
      else if (b==5) ret =  (ijk[2] == topo.idx_max_for_level(level,2));
      else throw -239;
    }
  else throw -345;
  return ret;
}



bool oOctree::next_ijk(int* ijk, int level) const
{
  const bool debug=false;
    
  int imax = topo.idx_max_for_level(level,0);
  int jmax = topo.idx_max_for_level(level,1);
  int kmax = topo.idx_max_for_level(level,2);
//   jmax = (Dim >= 2)?imax:0;
  kmax = (Dim == 3)?kmax:0;


  if(debug) cout<<"level="<<level<<" maxi X "<<imax<<" Y "<<jmax<<" Z "<<kmax<<endl;

  if (ijk[0] < imax)
    {
      ijk[0]++; return true;
    }
  else if (ijk[1] < jmax)
    {
      ijk[1]++; ijk[0]=0; return true;
    }
  else if (ijk[2] < kmax)
    {
      ijk[2]++; ijk[1]=0; ijk[0]=0; return true;
    }
  else return false;
}


void oOctree::coarsen_incompatible(oRefinementCriteria& refinement_criteria)
{
  const bool debug = false;
  if (debug) cout << "stating coarsen_incompatible" << endl;
  //we activate the most refined level
  //two lines below are now in the constructor
  //fill(begin(0),     end(LevelMax-1), 0);
  //fill(begin(LevelMax), end(LevelMax), 1);

  //coarsening stage
  for (int l = LevelMax -1; l >= 0; --l)
  {
    int ijk[3];
    fill(ijk, ijk+3, 0);
    do
    {
      //get the children
      iterator cell = cartesian2octree(ijk, l);
      iterator children = getChildren(cell, l);
      //check if all children are active.
      if (debug) std::cout << "NbChildren " << NbChildren << std::endl;

      iterator res = std::find_if(children, children+NbChildren, bind2nd(not_equal_to<int>(), 1));
      bool all_children_active = (res == children + NbChildren);

      if (debug) std::cout << "all_children_active: " << all_children_active << std::endl;

                        if (all_children_active)
            {
              if (!refinement_criteria(cell, l, ijk, children, children+NbChildren))  activate(cell, l, ijk);
            }
                }
                while ( next_ijk(ijk, l) );
    }
  if (debug) cout << "ending coarsen_incompatible" << endl;
}

void oOctree::refine_incompatible(oRefinementCriteria& refinement_criteria, int ldeb, int lfin)
{
  const bool debug = false;
  //refining stage
  if (ldeb >= lfin)
    {
      for (int l = ldeb; l >= lfin; --l)
        {
          if (debug) cout<<"l="<<l<<endl;
          int ijk[3];
          fill(ijk, ijk+3, 0);
          do
            {
              iterator cell = cartesian2octree(ijk, l);
              iterator children = getChildren(cell, l);
              if (*cell == 1) // we consider only active cells
                {
                  if(refinement_criteria(cell, l, ijk, children, children+NbChildren))  deactivate(cell, l, ijk);
                }
            }
          while ( next_ijk(ijk, l) );
        }

    }
  else
    {
      for (int l = ldeb; l <= lfin; ++l)
        {
          if (debug) cout<<"l="<<l<<endl;
          int ijk[3];
          fill(ijk, ijk+3, 0);
          do
            {
              iterator cell = cartesian2octree(ijk, l);
              iterator children = getChildren(cell, l);
              if (*cell == 1) // we consider only active cells
                {
                  if(refinement_criteria(cell, l, ijk, children, children+NbChildren))  deactivate(cell, l, ijk);
                }
            }
          while ( next_ijk(ijk, l) );
        }
    }
}

void oOctree::filter_cells(oFilterCriteria& filter_criteria)
{
  // loop over cells
  int ijk[3];
  for (int l = 0; l < LevelMax; ++l)
  {

    fill(ijk, ijk+3, 0);
    do
    {
      //get cell and his children
      iterator cell = cartesian2octree(ijk, l);
      iterator children = getChildren(cell, l);
      if (1 == *cell && filter_criteria(cell, l, ijk, children, children+NbChildren))
      {
//        iterator children = getChildren(cell, l);//not used, and dangerous (variable already exist...)
        deactivate(cell, l, ijk);
      }

    }
    while ( next_ijk(ijk, l) );
  }

  // deactivate the last level
  fill(ijk, ijk+3, 0);
  do
  {
    iterator cell = cartesian2octree(ijk, LevelMax);
    iterator children = nullptr;
    if (1 == *cell && filter_criteria(cell,LevelMax,ijk,children,
        children+NbChildren)) {
      *cell -= 1;
    }
  }
  while ( next_ijk(ijk, LevelMax) );

}


void oOctree::deactivateFinestLevel(){
  fill(begin(LevelMax), end(LevelMax), 0);
}

void oOctree::activateFinestLevel(){
  fill(begin(LevelMax), end(LevelMax), 1);
}


void oOctree::optimize(oRefinementCriteria& refinement_criteria, bool twoForOneRule)
{
  const bool debug = false;
  coarsen_incompatible(refinement_criteria);
  if (debug) cout << "in oOctree::optimize, starting refineForOneLevelDiffConstraint();" << endl;
  if(twoForOneRule) refineForOneLevelDiffConstraint();
}

void oOctree::refine(oRefinementCriteria& refinement_criteria)
{
  const bool debug = false;
  refine_incompatible(refinement_criteria, LevelMax-1, 0);
  if (debug) cout << "in oOctree::refine, starting refineForOneLevelDiffConstraint();" << endl;
  refineForOneLevelDiffConstraint();
}
void oOctree::refine(oRefinementCriteria& refinement_criteria, int ldeb, int lfin)
{
  refine_incompatible(refinement_criteria, ldeb, lfin);
  refineForOneLevelDiffConstraint();
}


void oOctree::refineForOneLevelDiffConstraint()
{
  const bool debug = false;
  const int on_what = 1;
  const int nb_neighbor = (Dim>2) ? 1 : -1;
  int neighbor;

  //refinement stage
  for (int l = LevelMax; l >= 0; --l)
  {
    int ijk[3];
    int ijk_up[3];
    fill(ijk, ijk+3, 0);
    do
    {
      iterator cell = cartesian2octree(ijk, l);
      if (*cell == 1) // we consider only active cells
          {
              bool flag;
              int count = 0;
              int ijkn[3];
          neighbor=0;
          while(topo.nextNeighbor(ijk, l, ijkn, on_what, flag, Dim, Periodicity, count, neighbor))
                  {
                            if (flag)
                                {
                                    cell_type* cneighbor = cartesian2octree(ijkn, l);
                                        for (int i = *cneighbor; i>2; i--)
                                        {
                                            int up = i-1;
                                                if (debug) cout<<"And ijkN="<<ijkn[0]<<" "<<ijkn[1]<<" "<<l<<endl;
                                                iterator ancestor = getAncestor(up, cneighbor, l);
                                                if (debug) cout << " ancestor activity is " << (int) *ancestor << endl;
                                                assert(*ancestor == 1);
                                                octree2cartesian (ancestor, l-up, ijk_up);
                                                deactivate(ancestor, l-up, ijk_up);
                                        }
                                }
                if (neighbor>nb_neighbor)
                {
                    neighbor=0;
                    count++;
                }
                else
                    ++neighbor;
                   }
            }
        }
        while ( next_ijk(ijk, l) );
  }

}





void oOctree::octree2cartesian  ( const cell_type* cell_octree, 
                                int level, int* ijk) const
{
  const bool debug = false;
  int Ioctree_base_nbchildren   [topo.MAX_DEPTH];
  int size_octree;

  const int Ioctree = cell_octree - begin(level);

  topo.decompose ( Ioctree, NbChildren, Ioctree_base_nbchildren, size_octree );

  if (debug)
  {
    std::cout << " coefs_octree  " << std::endl;
    std::copy(Ioctree_base_nbchildren, Ioctree_base_nbchildren+size_octree, std::ostream_iterator<int>(std::cout, " "));
    std::cout << endl;
  }

  std::fill(ijk, ijk+3, 0);

  int c = 1;
  for (int l=size_octree-1;l>=0;--l)
  {
    for (int d = 0; d < Dim; ++d)
    {
      ijk[d] += c * topo.base2_ijk[Ioctree_base_nbchildren[l]][d];
    }
    c*= 2;
  }

  if (debug)
  {
    cout << " ijk  " << endl;
    std::copy(ijk, ijk+3, std::ostream_iterator<int>(cout, " "));
    cout << endl;
  }
}


oOctree::cell_type* oOctree::cartesian2octree  ( const int* ijk, 
                                                 int level )const
{
  const bool debug = false;
  int Ioctree = ijk2_almost_octree[ijk[0]] +  2 * ijk2_almost_octree[ijk[1]];
  if (Dim == 3) Ioctree += 4 * ijk2_almost_octree[ijk[2]];
  if (debug)
    {
      cout << " in cartesian2octree ijk is " << endl;
      std::copy(ijk, ijk+3, std::ostream_iterator<int>(cout, " "));
      cout << " at level " << level << endl << " and Ioctree is " << Ioctree << endl;
    }
  return (oOctree::cell_type*)(begin(level) + Ioctree);
}


int oOctree::getNbActiveCells() const 
{
  const bool debug = false;
  if (debug) cout << " char " << sizeof(char) << " int " << sizeof(int) << endl;
  int nb = 0;
  int nb_tot =0;
  for (int l = 0; l <= LevelMax; ++l)
    {
      nb+= count(begin(l), end(l), 1);
      nb_tot += size(l);
    }
  if (debug) cout << " nb active " << nb << " nb tot " << nb_tot << endl;
  return nb;
}


  void oOctree::compress_descendants(const cell_type* cell , int level , std::set<unsigned int>& compressed) const
  {
    const bool debug = false;

    std::vector<cell_type> cells_local;
    
    iterator c = (iterator) cell;
    int size = 1;
    cells_local.push_back(*cell);
    for (int l = level; l < LevelMax; ++l)
      {
        size *= NbChildren;
        c = getChildren(c, l);
        cells_local.insert(cells_local.end(), c, c+size);
      }
    
    if (debug)
      {
        cout << " cells_local activity " << endl;
        std::copy(cells_local.begin(), cells_local.end(), std::ostream_iterator<int>(std::cout, " "));
        cout << endl;
      }
    oOctree::compress(cells_local.begin(), cells_local.end(), compressed);

  }

void oOctree::attach(oObserver& obs_)
{
  obs.push_back(&obs_);
}
void oOctree::detach(oObserver& obs_)
{
  obs.remove(&obs_);
}
void oOctree::notify_refinement(const cell_type* c, int l, const int* ijk)
{
  const bool debug = false;
  if (debug)
    {
      cout << " notify refinement for ijk : " << endl;
      std::copy(ijk, ijk+3, std::ostream_iterator<int>(cout, " ")); cout << " at level " << l << endl;;
    }

  std::list<oObserver*>::const_iterator it = obs.begin(), ite = obs.end();
  for ( ; it != ite; ++it)
    {
      (*it)->refine(*this, c, l, ijk);
    }
}
void oOctree::notify_derefinement(const cell_type* c, int l, const int* ijk)
{
  const bool debug = false;
  if (debug)
    {
      cout << " notify derefinement for ijk : " << endl;
      std::copy(ijk, ijk+3, std::ostream_iterator<int>(cout, " ")); cout << " at level " << l << endl;;
    }

  std::list<oObserver*>::const_iterator it = obs.begin(), ite = obs.end();
  for ( ; it != ite; ++it)
    {
      (*it)->derefine(*this, c, l, ijk);
    }
}



void oOctree::lookForNodeOnNeighbors2D(const cell_type* cell, int level, 
                                       const int* ijk, int* exist_on_edge) const
{
  const bool debug = false;
  fill(exist_on_edge,  exist_on_edge+4, 0);
  //loop over the edges
  int on_what = 1;
  int count = 0;
  int neighbor = 0;
  bool exist;
  int ijkn[3];
  while(topo.nextNeighbor(ijk, level, ijkn, on_what, exist, Dim, Periodicity, count, neighbor))
  {
    if (debug)
        {
        cout << " exist " << exist << endl;
        if (exist)
        {
                oOctree::cell_type* neighbor = cartesian2octree(ijkn, level);
                if (*neighbor == 0) cout << " neighbor activity is 0" << endl;
        }
        }

    if (exist && (*cartesian2octree(ijkn, level) == 0))
        {
            exist_on_edge[count] = 1;
        }
    count++;
  }
  if (debug)
  {
      cout << "action edge ";
      std::copy(exist_on_edge, exist_on_edge+4, std::ostream_iterator<double>(std::cout, " ")); cout << endl;
  }
}



void oOctree::lookForNodeOnNeighbors3D(const cell_type* cell, 
                                       int level, const int* ijk,
                                       int* exist_on_edge, int* exist_on_face) const
{
  const bool debug = false;
  if (debug)
  {
      cout << "lookForNodeOnNeighbors3D for ijk ";
      std::copy(ijk, ijk+3, std::ostream_iterator<int>(std::cout, " ")); cout << endl;
      cout << "on level "<<level<<" with cell status "<<(int)*cell<<endl;
  }

  fill(exist_on_edge,  exist_on_edge+12, 0);
  fill(exist_on_face,  exist_on_face+6, 0);
  //loop over the faces
  int on_what = 2;
  int count = 0;
  int neighbor = 0;
  bool exist;
  int ijkn[3];
  while(topo.nextNeighbor(ijk, level, ijkn, on_what, exist, Dim, Periodicity, count, neighbor))
  {
        if (exist && (*cartesian2octree(ijkn, level) == 0))
        {
            exist_on_face[count] = 1;
            if (debug)
            {
                cout << "ijkn ";
                std::copy(ijkn, ijkn+3, std::ostream_iterator<int>(std::cout, " ")); cout << endl;
                cout << "count "<<count<<" exist "<<exist<<" exist_on_face "<<exist_on_face[count]<<endl;
            }

            // this is a quick way to mark edges. It is used in loop on edge below to skip these already analysed edges
            for (int ed = 0; ed < 4; ++ed)
            {
                exist_on_edge[topo.face_edge_connect[count][ed]] = 1;
            }

        }
        count++;
  }
  //loop over the edges
  on_what = 1;
  count = 0;
  while(topo.nextNeighbor(ijk, level, ijkn, on_what, exist, Dim, Periodicity, count, neighbor))
    {
        if (debug)
        {
            cout << "ijkn ";
            std::copy(ijkn, ijkn+3, std::ostream_iterator<int>(std::cout, " ")); cout << endl;
            printf("count %d neighbor %d exist %d cartesian2octree (%ld) = %d",count,neighbor,exist,cartesian2octree(ijkn, level)-cells,*cartesian2octree(ijkn, level));
        }
        if (exist_on_edge[count]) neighbor=2;
        else if (exist && (*cartesian2octree(ijkn, level) == 0))
        {
            exist_on_edge[count] = 1;
            neighbor=2;
        }
        if (debug)
        {
            cout <<" exist_on_edge "<<exist_on_edge[count]<<endl;
        }
        if (neighbor>1)
        {
            neighbor=0;
            count++;
        }
        else
            ++neighbor;
    }
  if (debug) cout << "end lookForNodeOnNeighbors3D \n";
}

void oOctree::lookForNodeHangingOnNeighbors2D(const cell_type* cell, 
                                       int level, const int* ijk,
                                       int* hanging_on_edge) const
{
  const bool debug = false;
  if (debug)
  {
      cout << "lookForNodeHangingOnNeighbors2D for ijk ";
      std::copy(ijk, ijk+3, std::ostream_iterator<int>(std::cout, " ")); cout << endl;
      cout << "on level "<<level<<" with cell status "<<(int)*cell<<endl;
  }

  fill(hanging_on_edge,  hanging_on_edge+4, 1);
  //loop over the edges
  int on_what = 1;
  int count = 0;
  int neighbor = 0;
  bool exist;
  int ijkn[3];
  while(topo.nextNeighbor(ijk, level, ijkn, on_what, exist, Dim, Periodicity, count, neighbor))
  {
        if (exist )
        {
            if (*cartesian2octree(ijkn, level) > 0) ;
            else
                hanging_on_edge[count] = 0;
        }
        else
                hanging_on_edge[count] = 0;

        if (debug)
        {
            cout << "count "<<count<<" exist "<<exist<<" hanging_on_edge "<<hanging_on_edge[count]<<endl;
        }

        count++;
  }
  if (debug) cout << "end lookForNodeHangingOnNeighbors2D \n";
}

void oOctree::lookForNodeHangingOnNeighbors3D(const cell_type* cell, 
                                       int level, const int* ijk,
                                       int* hanging_on_edge, int* hanging_on_face) const
{
  const bool debug = false;
  if (debug)
  {
      cout << "lookForNodeHangingOnNeighbors3D for ijk ";
      std::copy(ijk, ijk+3, std::ostream_iterator<int>(std::cout, " ")); cout << endl;
      cout << "on level "<<level<<" with cell status "<<(int)*cell<<endl;
  }

  fill(hanging_on_edge,  hanging_on_edge+12, 3);
  fill(hanging_on_face,  hanging_on_face+6, 1);
  //loop over the faces
  int on_what = 2;
  int count = 0;
  int neighbor = 0;
  bool exist;
  int ijkn[3];
  while(topo.nextNeighbor(ijk, level, ijkn, on_what, exist, Dim, Periodicity, count, neighbor))
  {
        if (exist )
        {
            if (*cartesian2octree(ijkn, level) > 0) ;
            else
                hanging_on_face[count] = 0;
        }
        else
                hanging_on_face[count] = 0;

        if (debug)
        {
            cout << "count "<<count<<" exist "<<exist<<" hanging_on_face "<<hanging_on_face[count]<<endl;
        }

        count++;
  }
  //loop over the edges
  on_what = 1;
  count = 0;
  while(topo.nextNeighbor(ijk, level, ijkn, on_what, exist, Dim, Periodicity, count, neighbor))
    {
        if (exist)
        {
            if (*cartesian2octree(ijkn, level) > 0)  ;
            else
             hanging_on_edge[count]-=1;
        }
        else
        {
             hanging_on_edge[count]-=1;
        }
        if (neighbor>1)
        {
            if (debug)
            {
                cout <<"count "<<count<<" hanging_on_edge "<<hanging_on_edge[count]<<endl;
            }
            neighbor=0;
            count++;
        }
        else
            ++neighbor;
    }
  if (debug) cout << "end lookForNodeHangingOnNeighbors3D \n";
}

  void oOctree::modifyOctree(oObserverRecorder &obs){

  const bool debug=false;
  detach(obs);//we do not want to record the backward operations !

  oObserverRecorder::record_type::iterator itB=obs.begin();
  oObserverRecorder::record_type::iterator itE=obs.end();

    for(;itB!=itE;++itB) {
      bool state= itB->second;
      int level= (itB->first).first;
      std::vector<int> ijk=(itB->first).second;

      if(debug){
        cout<<"state="<<state<<" level="<<level<<" ijk="<<ijk[0]<<" "<<ijk[1]<<" "<<ijk[2]<<endl;
      }

      cell_type* cell=cartesian2octree(&ijk[0],level);
      if(state) activate(cell,level,&ijk[0]);
      else deactivate(cell,level,&ijk[0]);
      }

  obs.clear();//nothing to record now
  attach(obs);

  }

  void oOctree::activateComp  (oOctree::cell_type* cell, int level, const int* ijk)
  {
    iterator c = cell;
    int size = 1;
    *cell += 1;
    for (int l = level; l < LevelMax; ++l)
      {
        size *= NbChildren;
        c = getChildren(c, l);
        transform(c, c + size, c, bind2nd(std::plus<int>(), 1+l-level));
      }
    notify_derefinement(cell, level, ijk);
  }

  void oOctree::compareOctree(oObserverRecorder &obs){

  const bool debug=false;
  detach(obs);//we do not want to record the backward operations !

    //desactivate all the octree...
    fill(begin(0), end(LevelMax), 0);

  oObserverRecorder::record_type::iterator itB=obs.begin();
  oObserverRecorder::record_type::iterator itE=obs.end();

    for(;itB!=itE;++itB) {
      bool state= itB->second;
      int level= (itB->first).first;
      std::vector<int> ijk=(itB->first).second;

      if(debug){
        cout<<"state="<<state<<" level="<<level<<" ijk="<<ijk[0]<<" "<<ijk[1]<<" "<<ijk[2]<<endl;
      }

      cell_type* cell=cartesian2octree(&ijk[0],level);
      if(state) activateComp(cell,level,&ijk[0]);
      *cell-=1;

      }

  obs.clear();//nothing to record now
  attach(obs);

  }


  void oOctree::reset(){
    //most refined level is active
    fill(begin(0),     end(LevelMax-1), 0);
    fill(begin(LevelMax), end(LevelMax), 1);
  }

void oObserverRecorder::refine(const oOctree& octree, oOctree::const_iterator cell, int level, const int* ijk)
{
  const bool debug=false;
  if(debug) cout<<"record ref ijk="<<ijk[0]<<" "<<ijk[1]<<endl;
  record.push_back(make_pair(make_pair(level,std::vector<int>(ijk,ijk+3)), true));
}

void oObserverRecorder::derefine(const oOctree& octree, oOctree::const_iterator cell, int level, const int* ijk)
{
  const bool debug=false;
  if(debug) cout<<"record deref ijk="<<ijk[0]<<" "<<ijk[1]<<endl;
  record.push_back(make_pair(make_pair(level,std::vector<int>(ijk,ijk+3)), false));
}


void oObserverRecorder::getSons(const oOctree &octree, std::set<std::vector<int> > &sons){

  #if 1
//   std::set<std::vector<int> > sons;

  record_type::iterator it=begin();
  record_type::iterator ite=end();

  for(;it!=ite;++it){
     std::pair<int, std::vector<int> > splitPair=it->first;
     bool refined=it->second;

     if(refined){
       oOctree::cell_type* splitCell=octree.cartesian2octree( &splitPair.second[0], splitPair.first);

       oOctree::cell_type* children=octree.getChildren(splitCell,splitPair.first);
       for(int i=0;i<octree.getNbChildren();++i){
//          bornFromRefinement.insert(children++);
//         std::vector<int> ijk(3,0);
//         octree.octree2cartesian(children++,splitPair.first+1, &ijk[0]);
//         octree.octree2finest_cartesian(children++,splitPair.first+1, &ijk[0]);
//          sons.insert(make_pair(splitPair.first+1, ijk));
//           sons.insert(ijk);

          int ijk[3];
          fill(ijk, ijk+3, 0);
          octree.octree2finest_cartesian(children++,splitPair.first+1, ijk);

          cout<<"Insert Son "<<ijk[0]<<" "<<ijk[1]<<endl;;
           sons.insert(std::vector<int>(&ijk[0],&ijk[3]));

       }
     }
  }

  cout<<"Sons size="<<sons.size()<<endl;
  #endif


}


void oObserverRecorder::getGrandSons(const oOctree &octree, std::set<std::vector<int> > &gsons){

  #if 1
//   std::set<std::vector<int> > sons;

  record_type::iterator it=begin();
  record_type::iterator ite=end();

  for(;it!=ite;++it){
     std::pair<int, std::vector<int> > splitPair=it->first;
     bool refined=it->second;

     if(refined){
       oOctree::cell_type* splitCell=octree.cartesian2octree( &splitPair.second[0], splitPair.first);

       oOctree::cell_type* children=octree.getChildren(splitCell,splitPair.first);
       for(int i=0;i<octree.getNbChildren();++i){

         oOctree::cell_type* gchildren=octree.getChildren(children++,splitPair.first+1);

         for(int j=0;j<octree.getNbChildren();++j){

          int ijk[3];
          fill(ijk, ijk+3, 0);
          octree.octree2finest_cartesian(gchildren++,splitPair.first+2, ijk);

          cout<<"Insert Grand Son "<<ijk[0]<<" "<<ijk[1]<<endl;;
           gsons.insert(std::vector<int>(&ijk[0],&ijk[3]));
         }
       }
     }
  }

  cout<<"Grand Sons size="<<gsons.size()<<endl;
  #endif


}

//  bool oObserverRecorder::findSons(oOctree::cell_type *cell, int level){
//    
//    
//    return false;
//  }
//  
//  void oObserverRecorder::fillSons(){
//    
//   record_type::iterator it=record.begin();
//   record_type::iterator ite=record.end();
//   
//   for(it;it!=ite;++it){
//      pair<int, std::vector<int> > splitPair=it->first;
//      bool refined=it->second;
//      
//      if(refined){
//        oOctree::cell_type* splitCell=octree.cartesian2octree( &splitPair.second[0], splitPair.first);
//        
//        oOctree::cell_type* children=octree.getChildren(splitCell,splitPair.first);
//        for(int i=0;i<octree.getNbChildren();++i){
// //           recordSons.insert(children++);
//        }
//      }
//   }
//    
//    
//  }


void oOctree::locatePointOnFinestLevel(double x, double y, double z, int *ijk){

  int ijk0[3]={0,0,0};
  double binf[3]={0.,0.,0.};
  double bsup[3]={0.,0.,0.};
//  double step[3]={0.,0.,0.};

  mapping.getBox(ijk0,0,binf,bsup);
  const double *step=mapping.getStep(LevelMax);

//  int ijkF[3]={floor((x-binf[0])/step[0]), floor((y-binf[1])/step[1]), floor((z-binf[2])/step[2])};
  ijk[0]=floor((x-binf[0])/step[0]);
  ijk[1]=floor((y-binf[1])/step[1]);
  ijk[2]=floor((z-binf[2])/step[2]);
}

void oOctree::locatePointOnActiveCells(double x, double y, double z, int *ijk, int &level){

  int ijkF[3]={0,0,0};
  this->locatePointOnFinestLevel(x,y,z,ijkF);
  int currentLevel=LevelMax;
  oOctree::cell_type* ancestor= this->cartesian2octree(ijkF,currentLevel);

//  cout<<"currlev is "<<currentLevel<<endl;

  while(!(*ancestor == 1) && currentLevel>=1){
    ancestor=this->getAncestor(1,ancestor,currentLevel);
    currentLevel-=1;
  }

  if(currentLevel<1){
    ijk[0]=0;
    ijk[1]=0;
    ijk[2]=0;
    level=0;
  }else{
    this->octree2cartesian(ancestor,currentLevel,ijk);
    level=currentLevel;
  }

//  oOctree::cell_type* ancestor= this->cartesian2octree(ijkF,LevelMax);
//  oOctree::cell_type* me=ancestor;

//  if(*ancestor==1) cout<<"Already active at finest level\n";
//  else cout<<"Me not active\n";

//  int nbGenerations=0;
//  cout<<"bef\n";
//  while(!(*ancestor == 1)){
//    ++nbGenerations;
//    cout<<"up from "<<nbGenerations<<" generations\n";
//    ancestor=this->getAncestor(nbGenerations,me,LevelMax);
//    }
//  cout<<"after\n";
//    currentLevel-=nbGenerations;
//    this->octree2cartesian(ancestor,currentLevel,ijk);
//    level=currentLevel;


}


//inline void oOctree::getActiveLineage (const_iterator cell, int level, int level_max, list<pair<cell_type *, int> > &lineage, int generation){
//    cell_type *childrens=getChildren(cell, level);
//    for(int ic=0;ic<NbChildren;++ic){
//        cell_type *child=childrens;
//        if(*child) lineage.push_back(make_pair(child,level+1));
//        else if(level<= LevelMax) getActiveLineage (childrens, level+1, level_max, lineage, generation+1);
//        ++childrens;
//    }
//}












}//end namespace



