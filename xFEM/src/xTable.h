/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _XTABLE_H_
#define _XTABLE_H_

#include <map>
#include <vector>

class xTable
{
  typedef std::map<int,double> mt;
  mt m;
  std::vector<int> dimensions;
  
public :
  void clear() {dimensions.clear();m.clear();}
  void setdim1d() {dimensions.clear();dimensions.push_back(0);}
  void setdim2d(int j){dimensions.clear();dimensions.push_back(j);dimensions.push_back(0);}
  void setdim3d(int j,int k) {dimensions.clear();dimensions.push_back(j);dimensions.push_back(k);dimensions.push_back(0);}
  void setdim4d(int j,int k, int l) {dimensions.clear();dimensions.push_back(j);dimensions.push_back(k);dimensions.push_back(l);dimensions.push_back(0);}
  void setdimnd(int nb,int * tab) {dimensions.clear();for (int i=0;i<nb;++i) dimensions.push_back(tab[i]);}
  void setdimnd(const std::vector<int> &vec){dimensions=vec;}
  int getnbdim() const {return dimensions.size();}
  int getdim(unsigned int dim) const {assert(dim<dimensions.size());return dimensions[dim];}
  double get(int i) const
  {
    int s=dimensions.size();
    if ((s<1) && (i)) return 0.0;
    else
      {
      mt::const_iterator p=m.find(i);
      if (p!=m.end()) 
	return p->second; 
      else 
	return 0.0;
      }
  }

  double get(int i,int j) const
  {
    int s=dimensions.size();
    if (s<2)
      if (j) return 0.0;
      else return get(i);
    else
      return get(i+j*dimensions[0]);
  }

  double get(int i,int j,int k) const
  {
    int s=dimensions.size();
    if (s<3)
      if (k) return 0.0;
      else return get(i,j);
    else
      return get(i+dimensions[0]*j+dimensions[0]*dimensions[1]*k);
  }

  double get(int i,int j,int k,int l) const
  {
    int s=dimensions.size();
    if (s<4)
      if (l) return 0.0;
      else return get(i,j,k);
    else
      return get(i+dimensions[0]*j+dimensions[0]*dimensions[1]*k+dimensions[0]*dimensions[1]*dimensions[2]*l);
  }

  void set(int i,double v)
  { int s=dimensions.size();
    if ((s<1) && (i)) return ;
    else m[i]=v;
  }

  void set(int i,int j,double v)
  {
    int s=dimensions.size();
    if (s<2)
      if (j) return ;
      else set(i,v);
    else
      set(i+j*dimensions[0],v);
  }

  void set(int i,int j,int k,double v)
  {
    int s=dimensions.size();
    if (s<3)
      if (k) return ;
      else set(i,j,v);
    else
      set(i+j*dimensions[0]+k*dimensions[0]*dimensions[1],v);
  }

  void set(int i,int j,int k,int l,double v)
  {
    int s=dimensions.size();
    if (s<4)
      if (l) return ;
      else set(i,j,k,v);
    else
      set(i+j*dimensions[0]+k*dimensions[0]*dimensions[1]+l*dimensions[0]*dimensions[1]*dimensions[2],v);
  }
};
 
#endif //_XTABLE_H_
