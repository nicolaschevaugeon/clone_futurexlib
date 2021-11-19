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

using namespace std;

namespace xoctree {

oOctreeAnisotropic::oOctreeAnisotropic ( const oMapping& m_, int _LevelMax, int* per)
    : oOctree(m_,_LevelMax, per,false), motif(topo.getMotif()), anisotropic((motif[0]+motif[1]+motif[2])!=0), layer1D_0(topo.pow_base2[motif[0]]), layer2D_0(topo.pow_base2[motif[0]+motif[1]]), generation_size(topo.generation_size[Dim-1]) 
  {  

    const bool debug = false;

    // Pour bien avoir le motif de base le plus simple
    if(Dim==3) assert(motif[0]*motif[1]*motif[2]==0);
    if(Dim==2) assert(motif[0]*motif[1]==0);

    if (debug) cout<<"lx0_="<<motif[0]<<" ly0_="<<motif[1]<<" lz0_="<<motif[2]<<endl;

    construct(motif);

  }
  
  oOctreeAnisotropic::~oOctreeAnisotropic()
  = default;
  

  
void oOctreeAnisotropic::octree2cartesian  ( const cell_type* cell_octree, 
				int level, int* ijk) const
{
  const bool debug = false;
  int Ioctree_base_nbchildren   [topo.MAX_DEPTH];
  int size_octree;

  const int Ioctree = cell_octree - begin(level);
  if(debug) cout<<"Ioctree="<<Ioctree<<" niveau="<<level<<endl;

  std::fill(ijk, ijk+3, 0);
  int ijkLoc[3];
  std::fill(ijkLoc, ijkLoc+3, 0);

  //   Cas anisotrope :
  if(anisotropic){

//     Numerotation pour i puis j puis k croissants  
//     par analogie avec base2_ijk dans oTopo (sauf qu'ici, i,j et k ne varient pas tous)
//     layer1D_0= nb de cellules selon x au niveau 0
//     layer2D_0= nb de cellules selon (x,y) au niveau 0

    if(level==0){ 

    // k = Ioctree/taille d'une couche selon (x,y)
    // j = Le reste/taille d'une couche selon x
    // i = le reste
    
    ijk[2]= (int) (Ioctree) / layer2D_0;
    ijk[1]= (int) ((Ioctree )- ijk[2] * layer2D_0 ) / layer1D_0;
    ijk[0]= ((Ioctree +1)- ijk[2] * layer2D_0 - ijk[1] * layer1D_0) -1;
    }else{

      int generation=generation_size[level];//taille d'une generation d'enfants

      // Indice DANS la cellule du motif de base
      int Ibr=Ioctree%generation;
      // Indice DE la cellule du motif de base
      int Ib=((Ioctree-Ibr)/generation);
 
	if(debug){
// 	  cout<<"ijk="<<ijk[0]<<" "<<ijk[1]<<endl;
	  cout<<"Indice base="<<Ib<<" Indice base relatif="<<Ibr<<endl;
	  }

    ijk[2]= (int) (Ib) / layer2D_0;
    ijk[1]= (int) ((Ib )- ijk[2] * layer2D_0 ) / layer1D_0;
    ijk[0]= ((Ib +1)- ijk[2] * layer2D_0 - ijk[1] * layer1D_0) -1;

    transform(ijk, ijk+3, ijk, bind2nd(std::multiplies<int>(),topo.pow_base2[level]));

    if(debug){
	cout << " ijk Base :";
	std::copy(ijk, ijk+3, std::ostream_iterator<int>(cout, " ")); 
	cout << endl;
     }

      // Chaque cellule du motif de base contient un octree.
      // La cellule globale courante a un indice relatif Ibr dans le motif de base
      // On procede comme dans le cas isotrope pour calculer ijk, 
      // mais relativement a la cellule du motif de base
      topo.decompose ( Ibr, NbChildren, Ioctree_base_nbchildren, size_octree );

      if(debug){
	std::cout << " Relative coefs_octree :" ;
	std::copy(Ioctree_base_nbchildren, Ioctree_base_nbchildren+size_octree, std::ostream_iterator<int>(std::cout, " ")); 
	std::cout << endl;
      }

  
  std::fill(ijkLoc, ijkLoc+3, 0);

  int c = 1;
  for (int l=size_octree-1;l>=0;--l)
    {      
      for (int d = 0; d < Dim; ++d) 
	{
	  ijkLoc[d] += c * topo.base2_ijk[Ioctree_base_nbchildren[l]][d]; 
	}
      c*= 2;
    }

       if(debug){
	cout << " ijkLoc :";
	std::copy(ijkLoc, ijkLoc+3, std::ostream_iterator<int>(cout, " ")); 
	cout << endl;
      }

    transform(ijk, ijk+3, ijkLoc, ijk, std::plus<int>());

    }
}else{



  topo.decompose ( Ioctree, NbChildren, Ioctree_base_nbchildren, size_octree );

  if (debug) 
    {  
      std::cout << " coefs_octree  " << std::endl;
      std::copy(Ioctree_base_nbchildren, Ioctree_base_nbchildren+size_octree, std::ostream_iterator<int>(std::cout, " ")); 
      std::cout << endl;
    }


  int c = 1;
  for (int l=size_octree-1;l>=0;--l)
    {      
      for (int d = 0; d < Dim; ++d) 
	{
	  ijk[d] += c * topo.base2_ijk[Ioctree_base_nbchildren[l]][d]; 
	}
      c*= 2;
    }

  }
  
  if (debug) 
    {
      cout << " ijk  " << endl;
      std::copy(ijk, ijk+3, std::ostream_iterator<int>(cout, " ")); 
      cout << endl;
    }
}


oOctree::cell_type* oOctreeAnisotropic::cartesian2octree  ( const int* ijk, 
						 int level )const 
{
  const bool debug = false;

  int Ioctree;

  if(1){

    if(level<1){
    Ioctree=ijk[0]+ ijk[1]*layer1D_0 + ijk[2]*layer2D_0;
    }else{

    // ijk DANS la cellule du motif de base
//     int ijkBr[3]={ijk[0]%topo.pow_base2[level],
// 		  ijk[1]%topo.pow_base2[level],
// 		  ijk[2]%topo.pow_base2[level]};

	int ijkBr[3];
	std::transform(ijk,ijk+3,ijkBr,bind2nd(std::modulus<int>(), topo.pow_base2[level]));

    // ijk DE la cellule du motif de base
//       int ijkB[3]={(ijk[0]-ijkBr[0])/topo.pow_base2[level],
// 		   (ijk[1]-ijkBr[1])/topo.pow_base2[level],
// 		   (ijk[2]-ijkBr[2])/topo.pow_base2[level]};

      int ijkB[3];
      std::transform(ijk,ijk+3,ijkBr,ijkB,std::minus<int>());
      std::transform(ijkB,ijkB+3,ijkB,bind2nd(std::divides<int>(), topo.pow_base2[level]));

    // Ioctree de la cellule du motif de base (numerotation consistante avec octree2cartesian)
    // car numerotation selon i, puis j et enfin k
    int IoctreeB=ijkB[0]+ ijkB[1]*layer1D_0 + ijkB[2]*layer2D_0;

    // Ioctree relativement a la cellule du motif de base
    int IoctreeR=ijk2_almost_octree[ijkBr[0]] +  2 * ijk2_almost_octree[ijkBr[1]];

    int generation=generation_size[level];//taille d'une generation d'enfants

    Ioctree = IoctreeB*generation + IoctreeR;//Ioctree final

    if (Dim == 3) Ioctree += 4 * ijk2_almost_octree[ijkBr[2]];

#if 0
    cout << " c2o ijkBr : " << endl;
    std::copy(ijkBr, ijkBr+3, std::ostream_iterator<int>(cout, " "));
    cout << endl;

    cout << " c2o ijkB : " << endl;
    std::copy(ijkB, ijkB+3, std::ostream_iterator<int>(cout, " "));
    cout << endl;
#endif


    }
  }else{

  Ioctree = ijk2_almost_octree[ijk[0]] +  2 * ijk2_almost_octree[ijk[1]];
  if (Dim == 3) Ioctree += 4 * ijk2_almost_octree[ijk[2]];
  }

  if (debug) 
    {
      cout << " in cartesian2octree ijk is " << endl;
      std::copy(ijk, ijk+3, std::ostream_iterator<int>(cout, " "));
      cout << " at level " << level << endl << " and Ioctree is " << Ioctree << endl;
    }
  return (oOctree::cell_type*)(begin(level) + Ioctree);
}


//void oOctreeAnisotropic::locatePointOnFinestLevel(double x, double y, double z, int *ijk){

//  int ijk0[3]={0,0,0};
//  double binf[3]={0.,0.,0.};
//  double bsup[3]={0.,0.,0.};
////  double step[3]={0.,0.,0.};

//  mapping.getBox(ijk0,0,binf,bsup);
//  const double *step=mapping.getStep(LevelMax);

////  int ijkF[3]={floor((x-binf[0])/step[0]), floor((y-binf[1])/step[1]), floor((z-binf[2])/step[2])};
//  ijk[0]=floor((x-binf[0])/step[0]);
//  ijk[1]=floor((y-binf[1])/step[1]);
//  ijk[2]=floor((z-binf[2])/step[2]);
//}

//void oOctreeAnisotropic::locatePointOnActiveCells(double x, double y, double z, int *ijk, int &level){

//  int ijkF[3]={0,0,0};
//  this->locatePointOnFinestLevel(x,y,z,ijkF);
//  int currentLevel=LevelMax;
//  oOctree::cell_type* ancestor= this->cartesian2octree(ijkF,currentLevel);

//  cout<<"currleva is "<<currentLevel<<endl;

//  while(!(*ancestor == 1) || currentLevel<1){
//    ancestor=this->getAncestor(1,ancestor,currentLevel);
//    currentLevel-=1;
//  }

//  if(currentLevel<1){
//    ijk[0]=0;
//    ijk[1]=0;
//    ijk[2]=0;
//    level=-1;
//  }else{
//    this->octree2cartesian(ancestor,currentLevel,ijk);
//    level=currentLevel;
//  }

//}


}//end namespace



