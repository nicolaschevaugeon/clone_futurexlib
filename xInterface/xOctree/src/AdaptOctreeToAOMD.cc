/* 
    xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE & LICENSE files for conditions.
*/

#include "xMesh.h"
#include "xLevelSet.h"
#include "oOctree.h"
#include "AdaptOctreeToAOMD.h"
#include "oKeyManager.h"
#include "oExport.h"

#include "mTet.h"
#include "mHex.h"
#include "mFace.h"
#include <algorithm>

//main
using namespace xfem;
using namespace AOMD;

#define ALTERNATE 0
using namespace xoctree;


//xz/normale -y -> surf 1
//xz/normale +y -> surf 3
//xy/normale -z -> surf 5
//xy/normale +z -> surf 6
//yz/normale -x -> surf 4
//yz/normale +x -> surf 2





namespace xinterface {

namespace xoctree {

//   void updateActiveNodes(const oOctree &octree,oActiveNodes& active_nodes)
//   {
//     //octree.refineForOneLevelDiffConstraint();/*Pour eviter les oublis...*/
//     active_nodes.clearActiveNodes();
//     active_nodes.initializeMap(octree);
//   }



void ReadLsetFromITK2D(const string filename, std::vector<double> &lsVec, std::vector<int> &sizeLs, std::vector<double> &stepPx, std::vector<double> &origin){

    int sizeX,sizeY;
    double spaceX,spaceY,originX,originY;

    std::ifstream infile(filename.c_str());
    infile>>sizeX;
    infile>>sizeY;
    infile>>spaceX;
    infile>>spaceY;
    infile>>originX;
    infile>>originY;


    cout<<"ReadLevelSet From ITK 2D :\n";
    cout<<"size x,y :"<<sizeX<<" "<<sizeY<<endl;
    cout<<"space x,y :"<<spaceX<<" "<<spaceY<<endl;
    cout<<"origin x,y :"<<originX<<" "<<originY<<endl;

    lsVec.reserve(size_t(sizeX*sizeY));

    double value =0.;

    for(int i=0;i<sizeX*sizeY;++i){
        infile>>value;
        lsVec.push_back(value);
    }

    sizeLs.resize(2);
    sizeLs[0]=sizeX;
    sizeLs[1]=sizeY;
    stepPx[0]=spaceX;
    stepPx[1]=spaceY;
    origin[0]=originX;
    origin[1]=originY;



}


void ReadLsetFromITK3D(const string filename, std::vector<double> &lsVec, std::vector<int> &sizeLs, std::vector<double> &stepPx, std::vector<double> &origin, bool invertSide){

    int dim,sizeX,sizeY,sizeZ;
    double spaceX,spaceY,spaceZ,originX,originY,originZ;

    std::ifstream infile(filename.c_str());
    infile>>dim;
    infile>>sizeX;
    infile>>sizeY;
    infile>>sizeZ;
    infile>>spaceX;
    infile>>spaceY;
    infile>>spaceZ;
    infile>>originX;
    infile>>originY;
    infile>>originZ;


    cout<<"ReadLevelSet From ITK 3D :\n";
    cout<<"size x,y,z :"<<sizeX<<" "<<sizeY<<" "<<sizeZ<<endl;
    cout<<"space x,y,z :"<<spaceX<<" "<<spaceY<<" "<<spaceZ<<endl;
    cout<<"origin x,y,z :"<<originX<<" "<<originY<<" "<<originZ<<endl;

    lsVec.reserve(size_t(sizeX*sizeY*sizeZ));

    double value =0.;

    for(int i=0;i<sizeX*sizeY*sizeZ;++i){
        infile>>value;
        if(invertSide) value*=-1.;
        lsVec.push_back(value);
    }

    sizeLs.resize(3);
    sizeLs[0]=sizeX;
    sizeLs[1]=sizeY;
    sizeLs[2]=sizeZ;
    stepPx.resize(3);
    stepPx[0]=spaceX;
    stepPx[1]=spaceY;
    stepPx[2]=spaceZ;
    origin.resize(3);
    origin[0]=originX;
    origin[1]=originY;
    origin[2]=originZ;



}



}//end subnamespace

}//end namespace


// ---------------------------------------
