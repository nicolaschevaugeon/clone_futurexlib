/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <iostream>
#include <cmath>
#include "xExport.h"
#include "xVector.h"
#include "xTensor2.h"


using std::cout;
using std::endl;

namespace xfem {

static void PascalgetIndices(int iFct, int &n, int &i)
{
  int k = 0;
  int l = 0;
  while(k<=iFct)
    {
      l++;
      k +=l;
    };
  n = l - 1;
  i = l - k + iFct;
}

/*-----EDGE Beg------*/

xSplitEdge::xSplitEdge(int l, const Trellis_Util::mPoint &p1, const Trellis_Util::mPoint &p2):level(l)
{
  thisEdge = 0;
  nF = 0;
  t[0] = p1;t[1] = p2;
}

int xSplitEdge:: nbEdge (){

  double res = level;
  //res=pow(double(2),level-1);
  return (int) res;
}

bool xSplitEdge::nextEdge(Trellis_Util::mPoint &pt1, Trellis_Util::mPoint &pt2)
{
  const bool debug=false;
  if(nF++ == nbEdge()){cout<<"FALSE\n"; return false;}
  
  // pt1 = Trellis_Util::mPoint(-1.+(nF-1)*2./nbEdge(),0.,0.);
  //pt2 = Trellis_Util::mPoint(-1.+(nF)*2./nbEdge(),0.,0.);
  // previous is changed to the following by nicolas chevaugeon
  pt1 = Trellis_Util::mPoint(t[0]+(t[1]-t[0])*(nF-1)*(1./nbEdge()));
  pt2 = Trellis_Util::mPoint(t[0]+(t[1]-t[0])*(nF)*(1./nbEdge()));
  if(debug){
    cout<<"nF="<<nF<<" nb Edges="<<nbEdge()<<endl;
    cout<<"PT1 SubEdge "<<pt1(0)<<"  "<<pt1(1)<<"  "<<pt1(2)<<endl;
    cout<<"PT2 SubEdge "<<pt2(0)<<"  "<<pt2(1)<<"  "<<pt2(2)<<endl;
  }
  return true;
}


/*-----EDGE End------*/


xSplitTri::xSplitTri(int l, const Trellis_Util::mPoint &p1, const Trellis_Util::mPoint &p2, const Trellis_Util::mPoint &p3):level(l)
{
  thisTri = 0;
  iF = nF = 0;
  t[0] = p1;t[1] = p2;t[2] = p3;
}

int xSplitTri:: nbTri (){
  return (level)*(level);
}

bool xSplitTri::nextTri(Trellis_Util::mPoint &pt1, Trellis_Util::mPoint &pt2, Trellis_Util::mPoint &pt3)
{
  int i,n;
  if(nF++ == level*level)return false;
  PascalgetIndices(iF,i,n);
  if(i == n || thisTri == 1)
    {
      iF++;
      thisTri = 0;
    }
  else
    {
      thisTri++;
    }
  
  double x1 = (double)i/(double)(level);
  double y1 = (double)n/(double)(level);
  double x2 = (double)(i+1)/(double)(level);
  double y2 = (double)(n+1)/(double)(level);
  if(thisTri == 0)
    {
      pt1 = t[1] + ((t[0]-t[1]) * x1) + ((t[2]-t[0]) * y1);
      pt2 = t[1] + ((t[0]-t[1]) * x2) + ((t[2]-t[0]) * y1);
      pt3 = t[1] + ((t[0]-t[1]) * x2) + ((t[2]-t[0]) * y2);
    }
  else
    {
      pt1 = t[1] + ((t[0]-t[1]) * x1) + ((t[2]-t[0]) * y1);
      pt2 = t[1] + ((t[0]-t[1]) * x1) + ((t[2]-t[0]) * y2);
      pt3 = t[1] + ((t[0]-t[1]) * x2) + ((t[2]-t[0]) * y2);
    }
  if(thisTri %2){
    Trellis_Util::mPoint ptmp;
    ptmp = pt1;
    pt1 = pt2;
    pt2 = ptmp;
  }
  return true;
}



void recur_split_quad ( const Trellis_Util::mPoint &p0,
					 const Trellis_Util::mPoint &p1,
					 const Trellis_Util::mPoint &p2,
					 const Trellis_Util::mPoint &p3,
					 std::vector<Trellis_Util::mPoint> lst[4],
					 int depth, int max_depth)
{
  depth++;
  if(depth >= max_depth)
    {
      lst[0].push_back(p0);
      lst[1].push_back(p1);
      lst[2].push_back(p2);
      lst[3].push_back(p3);
    }
  else
    {
      Trellis_Util::mPoint pe0 = (p0+p1) * 0.5;
      Trellis_Util::mPoint pe1 = (p1+p2) * 0.5;
      Trellis_Util::mPoint pe2 = (p2+p3) * 0.5;
      Trellis_Util::mPoint pe3 = (p3+p0) * 0.5;
      Trellis_Util::mPoint pe4 = (pe3+pe1) * 0.5;
      recur_split_quad (p0 ,pe0,pe4,pe3,lst,depth,max_depth);
      recur_split_quad (pe0 ,p1,pe1,pe4,lst,depth,max_depth);
      recur_split_quad (pe3 ,pe4,pe2,p3,lst,depth,max_depth);
      recur_split_quad (pe4 ,pe1,p2,pe2,lst,depth,max_depth);
    }
}

xSplitQuad::xSplitQuad(int l, const Trellis_Util::mPoint &p1, const Trellis_Util::mPoint &p2, const Trellis_Util::mPoint &p3, const Trellis_Util::mPoint &p4):level(l)
{
  thisQuad = 0;
  recur_split_quad(p1,p2,p3,p4,v,0,l);
}

int xSplitQuad:: nbQuad (){
  return v[0].size();
}

bool xSplitQuad::nextQuad(Trellis_Util::mPoint &pt1, Trellis_Util::mPoint &pt2, Trellis_Util::mPoint &pt3, Trellis_Util::mPoint &pt4)
{
  if(thisQuad == nbQuad())
    {
      thisQuad = 0;
      return false;
    }
  pt1 = v[0][thisQuad];
  pt2 = v[1][thisQuad];
  pt3 = v[2][thisQuad];
  pt4 = v[3][thisQuad];
  thisQuad++;
  return true;
}

/*** SPLIT TETS *****/

void recur_split_tet ( const Trellis_Util::mPoint &p0,
		       const Trellis_Util::mPoint &p1,
		       const Trellis_Util::mPoint &p2,
		       const Trellis_Util::mPoint &p3, 
		       std::vector<Trellis_Util::mPoint> lst[4],
		       int depth, int max_depth)
{
  depth++;
  if(depth >= max_depth)
    {
      lst[0].push_back(p0);
      lst[1].push_back(p1);
      lst[2].push_back(p2);
      lst[3].push_back(p3);
    }
  else
    {
      Trellis_Util::mPoint pe0 = (p0+p1) * 0.5;
      Trellis_Util::mPoint pe1 = (p0+p2) * 0.5;
      Trellis_Util::mPoint pe2 = (p0+p3) * 0.5;
      Trellis_Util::mPoint pe3 = (p1+p2) * 0.5;
      Trellis_Util::mPoint pe4 = (p1+p3) * 0.5;
      Trellis_Util::mPoint pe5 = (p2+p3) * 0.5;
      recur_split_tet (p0 ,pe0,pe2,pe1,lst,depth,max_depth); 
      recur_split_tet (p1 ,pe0,pe3,pe4,lst,depth,max_depth); 
      recur_split_tet (p2 ,pe3,pe1,pe5,lst,depth,max_depth); 
      recur_split_tet (p3 ,pe2,pe4,pe5,lst,depth,max_depth); 
      recur_split_tet (pe3,pe5,pe2,pe4,lst,depth,max_depth); 
      recur_split_tet (pe3,pe2,pe0,pe4,lst,depth,max_depth); 
      recur_split_tet (pe2,pe5,pe3,pe1,lst,depth,max_depth); 
      recur_split_tet (pe0,pe2,pe3,pe1,lst,depth,max_depth); 
    }
}

xSplitTet::xSplitTet(int l, const Trellis_Util::mPoint &p1, const Trellis_Util::mPoint &p2, const Trellis_Util::mPoint &p3 , const Trellis_Util::mPoint &p4)
{
  thisTet = 0;
  recur_split_tet(p1,p2,p3,p4,v,0,l);
}

int xSplitTet:: nbTet (){
  return v[0].size();
}

bool xSplitTet::nextTet(Trellis_Util::mPoint &pt1, Trellis_Util::mPoint &pt2, Trellis_Util::mPoint &pt3, Trellis_Util::mPoint &pt4)
{
  if(thisTet == nbTet())
    {
      thisTet = 0;
      return false;
    }
  pt1 = v[0][thisTet];
  pt2 = v[1][thisTet];
  pt3 = v[2][thisTet];
  pt4 = v[3][thisTet];
  thisTet++;
  return true;
}

/*** SPLIT HEX *****/


void recur_split_hex ( const Trellis_Util::mPoint &p0,
											 const Trellis_Util::mPoint &p1,
											 const Trellis_Util::mPoint &p2,
											 const Trellis_Util::mPoint &p3,
											 const Trellis_Util::mPoint &p4,
											 const Trellis_Util::mPoint &p5,
											 const Trellis_Util::mPoint &p6,
											 const Trellis_Util::mPoint &p7,
											 std::vector<Trellis_Util::mPoint> lst[8],
											 int depth, int max_depth)
{
  depth++;
  if(depth >= max_depth)
    {
      lst[0].push_back(p0);
      lst[1].push_back(p1);
      lst[2].push_back(p2);
      lst[3].push_back(p3);
      lst[4].push_back(p4);
      lst[5].push_back(p5);
      lst[6].push_back(p6);
      lst[7].push_back(p7);
    }
  else
    {
      Trellis_Util::mPoint pe0 = (p0+p1) * 0.5;
      Trellis_Util::mPoint pe1 = (p1+p5) * 0.5;
      Trellis_Util::mPoint pe2 = (p4+p5) * 0.5;
      Trellis_Util::mPoint pe3 = (p4+p0) * 0.5;
      Trellis_Util::mPoint pe4 = (p1+p2) * 0.5;

      Trellis_Util::mPoint pe5 = (p2+p6) * 0.5;
      Trellis_Util::mPoint pe6 = (p5+p6) * 0.5;
      Trellis_Util::mPoint pe7 = (p6+p7) * 0.5;
      Trellis_Util::mPoint pe8 = (p3+p7) * 0.5;
      Trellis_Util::mPoint pe9 = (p3+p2) * 0.5;

      Trellis_Util::mPoint pe10 = (p4+p7) * 0.5;
      Trellis_Util::mPoint pe11 = (p0+p3) * 0.5;
      Trellis_Util::mPoint pe12 = (pe0+pe2) * 0.5;
      Trellis_Util::mPoint pe13 = (pe1+pe5) * 0.5;
      Trellis_Util::mPoint pe14 = (pe5+pe8) * 0.5;

      Trellis_Util::mPoint pe15 = (pe3+pe8) * 0.5;
      Trellis_Util::mPoint pe16 = (pe2+pe7) * 0.5;
      Trellis_Util::mPoint pe17 = (pe0+pe9) * 0.5;
      Trellis_Util::mPoint pe18 = (pe17+pe16) * 0.5;

      recur_split_hex (p0 ,pe0,pe17,pe11,pe3,pe12,pe18,pe15, lst,depth,max_depth);
      recur_split_hex (pe0 ,p1,pe4,pe17,pe12,pe1,pe13,pe18, lst,depth,max_depth);
      recur_split_hex (pe17 ,pe4,p2,pe9,pe18,pe13,pe5,pe14, lst,depth,max_depth);
      recur_split_hex (pe11 ,pe17,pe9,p3,pe15,pe18,pe14,pe8, lst,depth,max_depth);
      recur_split_hex (pe3 ,pe12,pe18,pe15,p4,pe2,pe16,pe10, lst,depth,max_depth);
      recur_split_hex (pe12 ,pe1,pe13,pe18,pe2,p5,pe6,pe16, lst,depth,max_depth);
      recur_split_hex (pe18 ,pe13,pe5,pe14,pe16,pe6,p6,pe7, lst,depth,max_depth);
      recur_split_hex (pe15 ,pe18,pe14,pe8,pe10,pe16,pe7,p7, lst,depth,max_depth);


    }
}



xSplitHex::xSplitHex(int l, const Trellis_Util::mPoint &p1, const Trellis_Util::mPoint &p2, const Trellis_Util::mPoint &p3, const Trellis_Util::mPoint &p4, const Trellis_Util::mPoint &p5, const Trellis_Util::mPoint &p6, const Trellis_Util::mPoint &p7, const Trellis_Util::mPoint &p8):level(l)
{
  thisHex = 0;
  recur_split_hex(p1,p2,p3,p4,p5,p6,p7,p8,v,0,l);
}

int xSplitHex:: nbHex (){
  return v[0].size();
}

bool xSplitHex::nextHex(Trellis_Util::mPoint &pt1, Trellis_Util::mPoint &pt2, Trellis_Util::mPoint &pt3, Trellis_Util::mPoint &pt4,Trellis_Util::mPoint &pt5, Trellis_Util::mPoint &pt6, Trellis_Util::mPoint &pt7, Trellis_Util::mPoint &pt8)
{
  if(thisHex == nbHex())
    {
      thisHex = 0;
      return false;
    }

  pt1 = v[0][thisHex];
  pt2 = v[1][thisHex];
  pt3 = v[2][thisHex];
  pt4 = v[3][thisHex];
  pt5 = v[4][thisHex];
  pt6 = v[5][thisHex];
  pt7 = v[6][thisHex];
  pt8 = v[7][thisHex];
  thisHex++;
  return true;
}

}
