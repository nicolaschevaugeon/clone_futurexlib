/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

  

#ifndef _ANN_INTERFACE_
#define _ANN_INTERFACE_

#ifdef HAVE_ANN
#include "ANN/ANN.h"
#endif

#include <vector>
#include "mVertex.h"
#include <iostream>

namespace xfem{
class xMesh;
}






/// xNearestNeighborInterface create a data structure using ANN, for an iterator range on mVertex *. a call to nearestvertex to then return the closest vertex o f the range to the input point.
template <class ITER>
class xNearestNeighborInterface{
 	public:
		xNearestNeighborInterface(ITER itb, ITER end){
#ifdef HAVE_ANN
			int maxPts =0;
			ITER it =itb;
			while (it!=end){++it;++maxPts;}
			varray.reserve(maxPts);
			it =itb;
			dataPts= annAllocPts(maxPts, 3);
			int i=0;
			while (it!=end){
				AOMD::mVertex * v= (AOMD::mVertex*) (*it);
				varray[i] = v;
				for (int k=0; k< 3; ++k)
					dataPts[i][k] = (v->point())(k);
				++it;
				++i;
			}
			anntree=new ANNkd_tree(dataPts, maxPts  ,  3);    
			return;
#else
			std::cerr << "Approximate Nearest Neighbor Library not defined. Please compile with the -HAVE_ANN=1 flag." << std::endl;
			throw;
#endif
		};

	 	AOMD::mVertex * nearestvertexto(const AOMD::mVertex &v, double & distance ){
	 		double xyz[3];
	 		for (int k=0; k< 3; ++k)
				xyz[k] = v.point()(k);
	 		return nearestvertexto(xyz, distance);
			//     ANNpoint q = annAllocPt(3);
			//     for (int k=0; k< 3; ++k)
			//       q[k] = v.point()(k);
			//     ANNidxArray  nn_idx = new ANNidx[1];   //
			//     ANNdistArray dists  = new ANNdist[1];;    //
			//     double       eps=0.0;
			//     anntree->annkSearch(q, 1, nn_idx, dists );
			//     distance = sqrt(dists [0]);
			//     //closestVertex = varray[ nn_idx[0]] ;
			//     int i = nn_idx[0];
			//     delete [] nn_idx;
			//     delete []  dists;
			//     return varray[i] ;
		}

	 	AOMD::mVertex * nearestvertexto(const double* xyz, double & distance ){
#ifdef HAVE_ANN
			ANNpoint q = annAllocPt(3);
			for (int k=0; k< 3; ++k)
				q[k] = xyz[k];
			ANNidxArray  nn_idx = new ANNidx[1];   //
			ANNdistArray dists  = new ANNdist[1];;    //
			//double       eps=0.0;
			anntree->annkSearch(q, 1, nn_idx, dists );
			distance = sqrt(dists [0]);
			//closestVertex = varray[ nn_idx[0]] ;
			int i = nn_idx[0];
			delete [] nn_idx;
			delete []  dists;
			annDeallocPt(q);
			return varray[i] ;
#else
			return nullptr;
#endif
		}

	 	std::vector<AOMD::mVertex*> nearestVerticesInsideRadius(const AOMD::mVertex &v, double radius_factor, std::vector<double>& distance ){
	 		double xyz[3];
	 		for (int k=0; k< 3; ++k)
				xyz[k] = v.point()(k);
	 		return nearestVerticesInsideRadius(xyz, radius_factor, distance);
		}

	 	std::vector<AOMD::mVertex*> nearestVerticesInsideRadius(const double* xyz, double radius_factor, std::vector<double>& distance ){
		std::vector<AOMD::mVertex*> vertices;
#ifdef HAVE_ANN
    ANNpoint q = annAllocPt(3);
    for (int k=0; k< 3; ++k)
      q[k] = xyz[k];
		// first, need to know the distance to the closest point
    ANNidxArray  nn_idx = new ANNidx[1];
    ANNdistArray dists  = new ANNdist[1];
    anntree->annkSearch(q, 1, nn_idx, dists );
    double radius = sqrt(dists[0]);
    delete [] nn_idx;
    delete []  dists;

		// now, computes nearest neighbors in radius
		double squared_distance = radius*radius*radius_factor*radius_factor;
    int number_of_points = anntree->annkFRSearch(q, squared_distance, 0);
    nn_idx = new ANNidx[number_of_points];
    dists  = new ANNdist[number_of_points];
    anntree->annkFRSearch(q, squared_distance, number_of_points, nn_idx, dists );
    for (int k=0;k<number_of_points;k++){
			distance.push_back(sqrt(dists[k]));
			vertices.push_back(varray[nn_idx[k]]);
		}
    delete [] nn_idx;
    delete []  dists;
		annDeallocPt(q);
#endif
    return vertices;
  }
  

	 	std::vector<AOMD::mVertex*> nearestVerticesInsideRadiusStrict(const AOMD::mVertex &v, double radius, std::vector<double>& distance ){
		std::vector<AOMD::mVertex*> vertices;
#ifdef HAVE_ANN
    ANNpoint q = annAllocPt(3);
    for (int k=0; k< 3; ++k)
      q[k] = v.point()(k);
		// computes nearest neighbors in radius
		double squared_distance = radius*radius;
    int number_of_points = anntree->annkFRSearch(q, squared_distance, 0);
    ANNidxArray nn_idx = new ANNidx[number_of_points];
    ANNdistArray dists  = new ANNdist[number_of_points];
    distance.reserve(number_of_points);
    vertices.reserve(number_of_points);
    anntree->annkFRSearch(q, squared_distance, number_of_points, nn_idx, dists );
    for (int k=0;k<number_of_points;k++){
			distance.push_back(sqrt(dists[k]));
			vertices.push_back(varray[nn_idx[k]]);
		}
    delete [] nn_idx;
    delete []  dists;
		annDeallocPt(q);
#endif
    return vertices;
  }
#ifdef HAVE_ANN
  ANNkd_tree *anntree; 
  ANNpointArray dataPts;
  std::vector<AOMD::mVertex *> varray;
#endif

  ~xNearestNeighborInterface(){
#ifdef HAVE_ANN
    annDeallocPts(dataPts);
    delete anntree;
    anntree=nullptr;
    annClose();
#endif
  }
};

/// Return the closest point to P on segment [AB]
inline void xClosestPointToEdge2d (const double *A, const double *B, const double *P, double *closestpointonedge){
  const double eps = 1.e-10;
  const double AB[2] = {B[0] - A[0], B[1]-A[1]};
  const double ABN = AB[0]*AB[0] + AB[1]*AB[1];
  const double AP[2] =  {P[0] - A[0], P[1]-A[1]};
  const double AP_AB =  AP[0]*AB[0] + AP[1]*AB[1];
  double u = (ABN< eps)? 0. : AP_AB/ABN;
  u = u< 0.?0:u;
  u = u> 1.?1.:u;
  closestpointonedge[0] = A[0] + u*AB[0];
  closestpointonedge[1] = A[1] + u*AB[1];
  return;
}

/// build a data structure from a 1d mesh on plane xy that is able to return, given a point, the closest point on this mesh. (not the closest node ... )  Brute Force is used : all edge are tested against one point. If it's to slow, a kd-tree for edge should be build first ... Future work
class xFindClosestTo1dMesh_BruteForce{
 public:
  xFindClosestTo1dMesh_BruteForce(const xfem::xMesh *_mesh);
  void nearestpointto(const double *P, double *Cp );
 private:
  const xfem::xMesh *mesh;
};


#endif
