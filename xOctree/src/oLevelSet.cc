/*
    octree is a subproject of  xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

#include "oLevelSet.h"
#include <iostream>
#include <numeric>

using namespace std;

namespace xoctree {


  oLevelSet::oLevelSet( const oMapping& m, int lmax)
    : mapping(m), level_max(lmax), dim(m.getDim()), topo(mapping.getTopo()) {}

  double oLevelSet::getLevelSet(const  int * ijk_node_fine) const
  { abort(); cerr << "not coded " << endl; return 0;}

  //  void oLevelSet::getLevelSetCurvatures(const  int * ijk, int level, std::vector<double> &curvatures) const
  //  {abort(); cerr << "getLevelSetCurvatures not coded " << endl; return;};

  double oLevelSet::getLevelSetCurvature(const  int * ijk_node, int level) const
  {abort(); cerr << "getLevelSetCurvature not coded " << endl; return 0;};

  void oLevelSet::getLevelSetCurvatures(const  int * ijk, int level, std::vector<double>& curvatures) const{

    int ijk_fine[3], ijk_n[3];
    topo.cartesian2finest_cartesian(ijk, level, ijk_fine, level_max);
    const int offset = topo.pow_base2[level_max-level];
    int count = 0;
    while(topo.next_ijk_node_on_element(ijk_n, 0, ijk_fine, offset, dim, count))
      {
        curvatures[count++] = getLevelSetCurvature(ijk_n, level_max);
      }

    return;
  }

  oLevelSetDiscrete::oLevelSetDiscrete(const oMapping& m, int lmax, double* ls_)
    : oLevelSet(m, lmax),   ls_ijk_fine(ls_), powp1(topo.pow_base2[level_max]+1) {}

  void oLevelSetDiscrete::getLevelSetValues (const int * ijk, int level,
                                             std::vector<double>& ls_vals) const
  {
    int ijk_fine[3], ijk_node[3];
    int c;
    topo.cartesian2finest_cartesian(ijk, level, ijk_fine, level_max);
    const int offset = topo.pow_base2[level_max-level];
    int count = 0;
    const int imaxp2 =  topo.index_max_for_levels[0][level_max]+2;
    const int imayp2 =  topo.index_max_for_levels[1][level_max]+2;
    while(topo.next_ijk_node_on_element(ijk_node, 0, ijk_fine, offset, dim, count))
      {
        c = ijk_node[0] + imaxp2 * ijk_node[1];
        if (dim == 3) c += imaxp2 * imayp2 * ijk_node[2];
        //        cout<<"C"<<c<<endl;
        ls_vals[count] = ls_ijk_fine[c];
        count++;
      }
    return;
  }
  
  oLevelSetAnalytical::oLevelSetAnalytical(const oMapping& m, int lmax, ls_xyz_t ls_)
    : oLevelSet(m, lmax),  ls_xyz(ls_) {}

  void oLevelSetAnalytical::getLevelSetValues (const int * ijk, int level,
                                               std::vector<double>& ls_vals) const
  {
    const bool debug = false;
    double inf[3], sup[3];
    mapping.getBox(ijk, level, inf, sup);
    int count = 0;
    if (dim >= 1) {
        ls_vals[count++] = ls_xyz(inf[0],inf[1],inf[2]);
        ls_vals[count++] = ls_xyz(sup[0],inf[1],inf[2]);
      }
    if (dim >= 2) {
        ls_vals[count++] = ls_xyz(sup[0],sup[1],inf[2]);
        ls_vals[count++] = ls_xyz(inf[0],sup[1],inf[2]);
      }
    if (dim == 3) {
        ls_vals[count++] = ls_xyz(inf[0],inf[1],sup[2]);
        ls_vals[count++] = ls_xyz(sup[0],inf[1],sup[2]);
        ls_vals[count++] = ls_xyz(sup[0],sup[1],sup[2]);
        ls_vals[count++] = ls_xyz(inf[0],sup[1],sup[2]);
      }
    if (debug)
      {
        cout << "export level set" << endl;
        std::copy(ls_vals.begin(), ls_vals.end(), std::ostream_iterator<double>(std::cout, " "));
        cout << endl;
      }
    return;
  }


  double oLevelSetDiscrete::getLevelSet(const int * ijk_node_fine) const
  {
    //    const int powp1 = topo.pow_base2[level_max]+1;
    //    int indice = ijk_node_fine[0] + powp1 * ijk_node_fine[1]+ powp1*powp1*ijk_node_fine[2];
    //    return ls_ijk_fine[indice];

    //    cout<<"--"<<ijk_node_fine[0]<<" "<<ijk_node_fine[1]<<" "<<getLevelSetIndex(ijk_node_fine)<<endl;
    return ls_ijk_fine[getLevelSetIndex(ijk_node_fine)];
  }
  
  double oLevelSetAnalytical::getLevelSet(const int * ijk_node_fine) const
  {
    double xyz[3]={0.,0.,0.};
    mapping.ijk2xyz(ijk_node_fine, level_max, xyz);
    return ls_xyz(xyz[0], xyz[1], xyz[2]);
  }



  bool oLevelSetDiscrete::is_node_on_boundary(int b, const int* ijk, int level) const
  {
    bool ret=false; // set to false to avoid -Wmaybe-uninitialized bellow
    if (dim == 2)
      {
        if (b == 0)    ret =  (ijk[1] == 0);
        else if (b==1) ret =  (ijk[0] == topo.idx_max_for_level(level,0)+1);
        else if (b==2) ret =  (ijk[1] == topo.idx_max_for_level(level,1)+1);
        else if (b==3) ret =  (ijk[0] == 0);
        else assert(0);
      }
    else if (dim == 3)
      {
        if (b == 0)    ret =  (ijk[1] == 0);
        else if (b==1) ret =  (ijk[0] == topo.idx_max_for_level(level,0)+1);
        else if (b==2) ret =  (ijk[1] == topo.idx_max_for_level(level,1)+1);
        else if (b==3) ret =  (ijk[0] == 0);
        else if (b==4) ret =  (ijk[2] == 0);
        else if (b==5) ret =  (ijk[2] == topo.idx_max_for_level(level,2)+1);
        else assert(0);
      }
    else assert(0);
    return ret;
  }

  int oLevelSetDiscrete::getLevelSetIndex(const int * ijk_node_fine) const
  {
    //    const int powp1 = topo.pow_base2[level_max]+1;
    return (ijk_node_fine[0] + powp1 * ijk_node_fine[1]+ powp1*powp1*ijk_node_fine[2]);
  }

  void oLevelSetDiscrete::getLevelSetIndices(const int * ijk, int level,
                                             std::vector<int>& ls_idx) const
  {
    int ijk_fine[3], ijk_node[3];
    int c;
    topo.cartesian2finest_cartesian(ijk, level, ijk_fine, level_max);
    const int offset = topo.pow_base2[level_max-level];
    int count = 0;
    const int imaxp2 =  topo.index_max_for_levels[0][level_max]+2;
    const int imayp2 =  topo.index_max_for_levels[1][level_max]+2;
    while(topo.next_ijk_node_on_element(ijk_node, 0, ijk_fine, offset, dim, count))
      {
        c = ijk_node[0] + imaxp2 * ijk_node[1];
        if (dim == 3) c += imaxp2 * imayp2 * ijk_node[2];
        ls_idx[count] = c;
        count++;
      }
    return;
  }


  //Hypothesis: ijk_node is on finest level
  double oLevelSetDiscrete::getLevelSetCurvature(const  int * ijk_node_fine, int level) const
  {
    const bool debug = false;
    if (debug)
      {
        cout << "in getLevelSetCurvature level is " << level << endl;
        cout << " ijk is ";
        std::copy(ijk_node_fine, ijk_node_fine+3, std::ostream_iterator<int>(cout, " "));
        cout << endl;
      }

    if(level != level_max) throw;

    int ijk_right_fine[3] = {ijk_node_fine[0]+1, ijk_node_fine[1], ijk_node_fine[2]};
    int ijk_left_fine[3] = {ijk_node_fine[0]-1, ijk_node_fine[1], ijk_node_fine[2]};
    int ijk_up_fine[3] = {ijk_node_fine[0], ijk_node_fine[1]+1, ijk_node_fine[2]};
    int ijk_down_fine[3] = {ijk_node_fine[0], ijk_node_fine[1]-1, ijk_node_fine[2]};
    int ijk_up_right_fine[3] = {ijk_node_fine[0]+1, ijk_node_fine[1]+1, ijk_node_fine[2]};
    int ijk_up_left_fine[3] = {ijk_node_fine[0]-1, ijk_node_fine[1]+1, ijk_node_fine[2]};
    int ijk_down_right_fine[3] = {ijk_node_fine[0]+1, ijk_node_fine[1]-1, ijk_node_fine[2]};
    int ijk_down_left_fine[3] = {ijk_node_fine[0]-1, ijk_node_fine[1]-1, ijk_node_fine[2]};


    //Check boundaries
    bool b0 = is_node_on_boundary(0,ijk_node_fine,level_max);
    bool b1 = is_node_on_boundary(1,ijk_node_fine,level_max);
    bool b2 = is_node_on_boundary(2,ijk_node_fine,level_max);
    bool b3 = is_node_on_boundary(3,ijk_node_fine,level_max);

    double curv;

    const double *step = mapping.getStep(level);
    double hx = step[0];
    double hy = step[1];
    double center = getLevelSet(ijk_node_fine);
    //    double right(0), left(0), up(0), down(0), up_right(0), up_left(0), down_right(0), down_left(0);

    //Do not evaluate values out of the domain: Correct ijk (We should rather use shifted operators in this case)
    if(b0){//Correct down's
        ijk_down_fine[1] += 1;
        ijk_down_right_fine[1] += 1;
        ijk_down_left_fine[1] += 1;
      }
    if(b1){//Correct right's
        ijk_right_fine[0] -= 1;
        ijk_up_right_fine[0] -= 1;
        ijk_down_right_fine[0] -= 1;
      }
    if(b2){//Correct up's
        ijk_up_fine[1] -= 1;
        ijk_up_right_fine[1] -= 1;
        ijk_up_left_fine[1] -= 1;
      }
    if(b3){//Correct left's
        ijk_left_fine[0] += 1;
        ijk_up_left_fine[0] += 1;
        ijk_down_left_fine[0] += 1;
      }


    double right  = getLevelSet(ijk_right_fine);
    double left   = getLevelSet(ijk_left_fine);
    double up     = getLevelSet(ijk_up_fine);
    double down   = getLevelSet(ijk_down_fine);
    double up_right = getLevelSet(ijk_up_right_fine);
    double up_left = getLevelSet(ijk_up_left_fine);
    double down_right = getLevelSet(ijk_down_right_fine);
    double down_left = getLevelSet(ijk_down_left_fine);


    //Centred operators: will give false results on the bnd. We should use shifted operators in this case...
    double phi_x = (right - left)/(2.*hx);
    double phi_y = (up - down)/(2.*hy);
    double phi_xx = (right-2.*center+left)/(hx*hx);
    double phi_yy = (up-2.*center+down)/(hy*hy);
    double phi_xy = (up_right - down_right  - (up_left - down_left) )/(4.*hx*hy);


    if (dim == 2)
      {
        double normGrad2 = phi_x*phi_x+phi_y*phi_y;

        if(normGrad2<1.e-12) curv = 0.;//Si pb symetrique/noeud -->cas pathologique avec norme du grad infinie...
        else{
            curv = (phi_xx*phi_y*phi_y-2.*phi_x*phi_y*phi_xy+
                    phi_yy*phi_x*phi_x)/pow(normGrad2,1.5);
          }
      }
    else{// DIM = 3
        double hz = step[2];

        int ijk_front_fine[3] = {ijk_node_fine[0], ijk_node_fine[1], ijk_node_fine[2]+1};
        int ijk_back_fine[3] = {ijk_node_fine[0], ijk_node_fine[1], ijk_node_fine[2]-1};

        int ijk_front_right_fine[3] = {ijk_node_fine[0]+1, ijk_node_fine[1], ijk_node_fine[2]+1};
        int ijk_back_right_fine[3] = {ijk_node_fine[0]+1, ijk_node_fine[1], ijk_node_fine[2]-1};
        int ijk_front_left_fine[3] = {ijk_node_fine[0]-1, ijk_node_fine[1], ijk_node_fine[2]+1};
        int ijk_back_left_fine[3] = {ijk_node_fine[0]-1, ijk_node_fine[1], ijk_node_fine[2]-1};

        int ijk_up_front_fine[3] = {ijk_node_fine[0], ijk_node_fine[1]+1, ijk_node_fine[2]+1};
        int ijk_up_back_fine[3] = {ijk_node_fine[0], ijk_node_fine[1]+1, ijk_node_fine[2]-1};
        int ijk_down_front_fine[3] = {ijk_node_fine[0], ijk_node_fine[1]-1, ijk_node_fine[2]+1};
        int ijk_down_back_fine[3] = {ijk_node_fine[0], ijk_node_fine[1]-1, ijk_node_fine[2]-1};

        //Check boundaries
        bool b4 = is_node_on_boundary(4,ijk_node_fine,level_max);
        bool b5 = is_node_on_boundary(5,ijk_node_fine,level_max);

        //        double front(0), back(0), front_right(0), back_right(0), front_left(0), back_left(0), up_front(0), up_back(0), down_front(0), down_back(0);


        //Do not evaluate values out of the domain: Correct ijk in this case (We should rather use shifted operators in this case)
        if(b0){//Correct down's
            ijk_down_front_fine[1] += 1;
            ijk_down_back_fine[1] += 1;
          }
        if(b1){//Correct right's
            ijk_front_right_fine[0] -= 1;
            ijk_back_right_fine[0] -= 1;
          }
        if(b2){//Correct up's
            ijk_up_front_fine[1] -= 1;
            ijk_up_back_fine[1] -= 1;
          }
        if(b3){//Correct left's
            ijk_front_left_fine[0] += 1;
            ijk_back_left_fine[0] += 1;
          }
        if(b4){//Correct back's
            ijk_back_fine[2] += 1;
            ijk_back_right_fine[2] += 1;
            ijk_back_left_fine[2] += 1;
            ijk_up_back_fine[2] += 1;
            ijk_down_back_fine[2] += 1;
          }
        if(b5){//Correct front's
            ijk_front_fine[2] -= 1;
            ijk_front_right_fine[2] -= 1;
            ijk_front_left_fine[2] -= 1;
            ijk_up_front_fine[2] -= 1;
            ijk_down_front_fine[2] -= 1;
          }




        double front = getLevelSet(ijk_front_fine);
        double back = getLevelSet(ijk_back_fine);

        //Centred operators: will give false results on the bnd. We should use shifted operators in this case...
        double phi_z = (front - back)/(2.*hz);
        double phi_zz = (front-2.*center+back)/(hz*hz);


        double front_right = getLevelSet(ijk_front_right_fine);
        double back_right = getLevelSet(ijk_back_right_fine);
        double front_left = getLevelSet(ijk_front_left_fine);
        double back_left = getLevelSet(ijk_back_left_fine);

        double up_front = getLevelSet(ijk_up_front_fine);
        double up_back = getLevelSet(ijk_up_back_fine);
        double down_front = getLevelSet(ijk_down_front_fine);
        double down_back = getLevelSet(ijk_down_back_fine);


        //Centred operators: will give false results on the bnd. We should use shifted operators in this case...
        double phi_xz = (front_right - back_right - front_left + back_left )/(4.*hx*hz);
        double phi_yz = (up_front - up_back - down_front + down_back )/(4.*hy*hz);

        double normGrad2 = phi_x*phi_x + phi_y*phi_y + phi_z*phi_z;

        if(normGrad2<1.e-12) curv = 0.;//Si pb symetrique/noeud -->cas pathologique avec norme du grad infinie...
        else{
            curv = (phi_xx * (phi_y*phi_y + phi_z*phi_z) + phi_yy * (phi_x*phi_x + phi_z*phi_z)
                    + phi_zz * (phi_x*phi_x + phi_y*phi_y)
                    -2.* (phi_xy*phi_x*phi_y + phi_xz*phi_x*phi_z + phi_yz*phi_y*phi_z)
                    )/pow(normGrad2,1.5);
          }


      }

    if (debug)
      {
        cout << " in AMR curvature is " << curv << " at node " << endl;
        copy(ijk_node_fine, ijk_node_fine+3, std::ostream_iterator<double>(cout, " "));
        cout << endl;
      }


    return curv;

  }






  

  double oLevelSetAnalytical::getLevelSetCurvature(const  int * ijk_node, int level) const
  {
    const bool debug = false;
    if (debug)
      {
        cout << "in getLevelSetCurvature level is " << level << endl;
        cout << " ijk is ";
        std::copy(ijk_node, ijk_node+3, std::ostream_iterator<int>(cout, " "));
        cout << endl;
      }
    double xyz[3];
    mapping.ijk2xyz(ijk_node, level, xyz);
    double curv;

    const double *step = mapping.getStep(level);
    double hx = step[0];
    double hy = step[1];
    double center = ls_xyz(xyz[0], xyz[1], xyz[2]);
    double right  = ls_xyz(xyz[0]+hx, xyz[1], xyz[2]);
    double left   = ls_xyz(xyz[0]-hx, xyz[1], xyz[2]);
    double up     = ls_xyz(xyz[0], xyz[1]+hy, xyz[2]);
    double down   = ls_xyz(xyz[0], xyz[1]-hy, xyz[2]);
    double up_right = ls_xyz(xyz[0]+hx, xyz[1]+hy, xyz[2]);
    double up_left = ls_xyz(xyz[0]-hx, xyz[1]+hy, xyz[2]);
    double down_right = ls_xyz(xyz[0]+hx, xyz[1]-hy, xyz[2]);
    double down_left = ls_xyz(xyz[0]-hx, xyz[1]-hy, xyz[2]);

    double phi_x = (right - left)/(2.*hx);
    double phi_y = (up - down)/(2.*hy);
    double phi_xx = (right-2.*center+left)/(hx*hx);
    double phi_yy = (up-2.*center+down)/(hy*hy);
    double phi_xy = (up_right - down_right  - (up_left - down_left) )/(4.*hx*hy);


    if (dim == 2)
      {
        double normGrad2 = phi_x*phi_x+phi_y*phi_y;
        //        cout<<"denom "<<normGrad2<<endl;

        if(normGrad2<1.e-12) curv = 0.;//Si pb symetrique/noeud -->cas pathologique avec norme du grad infinie...
        else{
            curv = (phi_xx*phi_y*phi_y-2.*phi_x*phi_y*phi_xy+
                    phi_yy*phi_x*phi_x)/pow(normGrad2,1.5);
          }

        //        cout<<curv<<endl;

      }
    // 		else assert(0); //to be coded for 3D
    else{// DIM = 3
        //        curv=0.;
        double hz = step[2];
        double front = ls_xyz(xyz[0], xyz[1], xyz[2]+hz);
        double back = ls_xyz(xyz[0], xyz[1], xyz[2]-hz);
        double phi_z = (front - back)/(2.*hz);
        double phi_zz = (front-2.*center+back)/(hz*hz);


        double front_right = ls_xyz(xyz[0]+hx, xyz[1], xyz[2]+hz);
        double back_right = ls_xyz(xyz[0]+hx, xyz[1], xyz[2]-hz);
        double front_left = ls_xyz(xyz[0]-hx, xyz[1], xyz[2]+hz);
        double back_left = ls_xyz(xyz[0]-hx, xyz[1], xyz[2]-hz);

        double up_front = ls_xyz(xyz[0], xyz[1]+hy, xyz[2]+hz);
        double up_back = ls_xyz(xyz[0], xyz[1]+hy, xyz[2]-hz);
        double down_front = ls_xyz(xyz[0], xyz[1]-hy, xyz[2]+hz);
        double down_back = ls_xyz(xyz[0], xyz[1]-hy, xyz[2]-hz);



        double phi_xz = (front_right - back_right - front_left + back_left )/(4.*hx*hz);
        double phi_yz = (up_front - up_back - down_front + down_back )/(4.*hy*hz);

        double normGrad2 = phi_x*phi_x + phi_y*phi_y + phi_z*phi_z;

        if(normGrad2<1.e-12) curv = 0.;//Si pb symetrique/noeud -->cas pathologique avec norme du grad infinie...
        else{
            curv = (phi_xx * (phi_y*phi_y + phi_z*phi_z) + phi_yy * (phi_x*phi_x + phi_z*phi_z)
                    + phi_zz * (phi_x*phi_x + phi_y*phi_y)
                    -2.* (phi_xy*phi_x*phi_y + phi_xz*phi_x*phi_z + phi_yz*phi_y*phi_z)
                    )/pow(normGrad2,1.5);
          }


      }
    
    if (debug)
      {
        cout << " in AMR curvature is " << curv << " at point " << endl;
        copy(xyz, xyz+3, std::ostream_iterator<double>(cout, " "));
        cout << endl;
      }


    return curv;

  }


  void oLevelSetAnalytical::convertToDiscrete(double* ls_) const{
    int imax = topo.idx_max_for_level(level_max,0)+1;
    int jmax = topo.idx_max_for_level(level_max,1)+1;
    int kmax = topo.idx_max_for_level(level_max,2)+1;
    if(dim == 2) kmax = 0;

    //    cout<<imax<<" "<<jmax<<" "<<kmax<<endl;

    const int powp2 = topo.idx_max_for_level(level_max,0)+2;
    double xyz[3]={0.,0.,0.};
    int ijk_node_fine[3];

    for(int k = 0 ; k < kmax+1 ; ++k){
        for(int j = 0 ; j < jmax+1 ; ++j){
            for(int i = 0 ; i < imax+1 ; ++i){
                ijk_node_fine[0] = i;
                ijk_node_fine[1] = j;
                ijk_node_fine[2] = k;

                //Si le format change pout la levelset discrete, on est dans la mouise...
                int idx = (ijk_node_fine[0] + powp2 * ijk_node_fine[1]+ powp2*powp2*ijk_node_fine[2]);
                //                cout<<i<<" "<<j<<" "<<k<<" idx "<<idx<<endl;
                mapping.ijk2xyz(ijk_node_fine, level_max, xyz);
                ls_[idx] = ls_xyz(xyz[0], xyz[1], xyz[2]);
              }
          }
      }
    return;
  }

  //  void oLevelSetAnalytical::getLevelSetCurvatures(const  int * ijk, int level, std::vector<double>& curvatures) const{

  //    int ijk_fine[3], ijk_n[3];
  //    topo.cartesian2finest_cartesian(ijk, level, ijk_fine, level_max);
  //    const int offset = topo.pow_base2[level_max-level];
  //    int count = 0;
  //    while(topo.next_ijk_node_on_element(ijk_n, 0, ijk_fine, offset, dim, count))
  //      {
  //        curvatures[count++] = getLevelSetCurvature(ijk_n, level_max);
  ////        double val = getLevelSetCurvature(ijk_n, level_max);
  ////        curvatures[count++] =val;
  ////        if(std::isnan(val)) {
  ////            cout<<"ijk n "<<ijk_n[0]<<" "<<ijk_n[1]<<" "<<val<<endl;
  ////            throw;
  ////          }
  //      }

  //    return;
  //  }

  oLevelSetOnActive::oLevelSetOnActive(const oOctree& o, const oField& f)
    : oLevelSet(o.getMapping(), o.getLevelMax()), octree(o), field(f), key_manager(f.getKeyManager())  {}

  void oLevelSetOnActive::getLevelSetValues (const int * ijk, int level,
                                             std::vector<double>& lsv) const
  {
    int ijk_fine[3], ijk_node[3];
    topo.cartesian2finest_cartesian(ijk, level, ijk_fine, level_max);
    const int offset = topo.pow_base2[level_max-level];
    int count = 0;
    while(topo.next_ijk_node_on_element(ijk_node, 0, ijk_fine, offset, dim, count))
      {
        const oKey* key = key_manager.find(ijk_node);
        assert(key);
        lsv[count] = field.getVal(key);
        count++;
      }
  }
  
  double oLevelSetOnActive::getLevelSet(const  int * ijk_node_fine) const
  {
    const oKey* key = key_manager.find(ijk_node_fine);
    if (!key) {
        cerr << " no value associated to location fine " <<  ijk_node_fine[0]
             << " " << ijk_node_fine[1] << " " << ijk_node_fine[2] << endl; abort();
      }
    return field.getVal(key);
  }

  //void oLevelSetOnActive::setLevelSet(const  int * ijk_node_fine, const double& val)
  //{
  //  keys.setVal( oKey(ijk_node_fine), val );
  //}
  
  oFilterCriteriaOnSign::oFilterCriteriaOnSign(const oMapping& m, const oLevelSet& ls_, bool revert_)
    : mapping(m), ls(ls_), dim(m.getDim()), revert(revert_)  {}
  bool oFilterCriteriaOnSign::operator()(oOctree::const_iterator cell, int level,
                                         const int * ijk,
                                         oOctree::const_iterator children_beg,
                                         oOctree::const_iterator children_end) const
  {
    std::vector<double> ls_vals(mapping.getTopo().pow_base2[dim]);
    ls.getLevelSetValues(ijk, level, ls_vals);

    if(!revert){
        double minValue = *std::min_element(ls_vals.begin(), ls_vals.end());
        return (minValue > 0.);
      }else{
        double maxValue = *std::max_element(ls_vals.begin(), ls_vals.end());
        return (maxValue < 0.);
      }
  }


  oFilterCriteriaCrossed::oFilterCriteriaCrossed(const oMapping& m, const oLevelSet& ls_, bool revert_)
    : mapping(m), ls(ls_), dim(m.getDim()), revert(revert_)  {}
  bool oFilterCriteriaCrossed::operator()(oOctree::const_iterator cell, int level,
                                          const int * ijk,
                                          oOctree::const_iterator children_beg,
                                          oOctree::const_iterator children_end) const
  {
    std::vector<double> ls_vals(mapping.getTopo().pow_base2[dim]);
    ls.getLevelSetValues(ijk, level, ls_vals);

    double minValue = *std::min_element(ls_vals.begin(), ls_vals.end());
    double maxValue = *std::max_element(ls_vals.begin(), ls_vals.end());

    if(revert) return (minValue*maxValue<=0.);
    else return (minValue*maxValue>0.);
  }
  
  oFilterCriteriaOnVolumeFraction::oFilterCriteriaOnVolumeFraction(const oMapping& m, const oLevelSet& ls_)
    : mapping(m), ls(ls_), dim(m.getDim())  {}
  bool oFilterCriteriaOnVolumeFraction::operator()(oOctree::const_iterator cell, int level,
                                                   const int * ijk,
                                                   oOctree::const_iterator children_beg,
                                                   oOctree::const_iterator children_end) const
  {
    std::vector<double> ls_vals(mapping.getTopo().pow_base2[dim]);
    ls.getLevelSetValues(ijk, level, ls_vals);
    double volumeFraction = 0.0;
    const unsigned int n = 4;
    if (dim == 2)
      {
        int NX = mapping.getTopo().pow_base2[n];
        int NY = mapping.getTopo().pow_base2[n];

        std::vector<double> levelSetValue((NX+1)*(NY+1), 0.0);
        int index = 0;
        double rate_x = 1.0/NX;
        double rate_y = 1.0/NY;
        for (int xx = 0; xx <= NX; xx++)
          for (int yy = 0; yy <= NY; yy++)
            {
              double ax = xx * rate_x;
              double ay = yy * rate_y;
              double value1 = ls_vals[0]*(1-ay) + ls_vals[3]*ay;
              double value2 = ls_vals[1]*(1-ay) + ls_vals[2]*ay;
              levelSetValue[index++] = value1*(1-ax) + value2*ax;
            }
        for (int nx = 1; nx <= NX; nx++)
          for (int ny = 1; ny <= NY; ny++)
            {
              double v1 = levelSetValue[nx-1 + (ny-1)*NX];
              double v2 = levelSetValue[nx + (ny-1)*NX];
              double v3 = levelSetValue[nx-1 + ny*NX];
              double v4 = levelSetValue[nx + ny*NX];
              if (v1 < 0 && v2 < 0 && v3 < 0 && v4 < 0) volumeFraction += 1;
              else if (v1 * v2 * v3 * v4 < 0) volumeFraction += 0.5;
            }
        
        volumeFraction /= NX * NY;
        //       std::cout << "volume fraction is " << volumeFraction << endl;
        
      }
    else if (dim == 3)
      {
        int NX = mapping.getTopo().pow_base2[n];
        int NY = mapping.getTopo().pow_base2[n];
        int NZ = mapping.getTopo().pow_base2[n];
        std::vector<double> levelSetValue((NX+1)*(NY+1)*(NZ+1), 0.0);
        int index = 0;
        double rate_x = 1.0/NX;
        double rate_y = 1.0/NY;
        double rate_z = 1.0/NZ;
        for (int xx = 0; xx <= NX; xx++)
          for (int yy = 0; yy <= NY; yy++)
            for (int zz = 0; zz <= NY; zz++)
              {
                double ax = xx * rate_x;
                double ay = yy * rate_y;
                double az = zz * rate_z;

                double valuez1 = ls_vals[0]*(1-az) + ls_vals[4]*az;
                double valuez2 = ls_vals[1]*(1-az) + ls_vals[5]*az;
                double valuez3 = ls_vals[2]*(1-az) + ls_vals[6]*az;
                double valuez4 = ls_vals[3]*(1-az) + ls_vals[7]*az;

                double valuey1 = valuez1*(1-ay) + valuez4*ay;
                double valuey2 = valuez2*(1-ay) + valuez3*ay;

                levelSetValue[index++] = valuey1*(1-ax) + valuey2*ax;
              }

        for (int nx = 1; nx <= NX; nx++)
          for (int ny = 1; ny <= NY; ny++)
            for (int nz = 1; nz <= NY; nz++)
              {
                double v1 = levelSetValue[nx-1 + (ny-1)*NX + (nz-1)*NX*NY];
                double v2 = levelSetValue[nx + (ny-1)*NX + (nz-1)*NX*NY];
                double v3 = levelSetValue[nx-1 + ny*NX + (nz-1)*NX*NY];
                double v4 = levelSetValue[nx + ny*NX + (nz-1)*NX*NY];

                double v5 = levelSetValue[nx-1 + (ny-1)*NX + nz*NX*NY];
                double v6 = levelSetValue[nx + (ny-1)*NX + nz*NX*NY];
                double v7 = levelSetValue[nx-1 + ny*NX + nz*NX*NY];
                double v8 = levelSetValue[nx + ny*NX + nz*NX*NY];

                if (v1 < 0 && v2 < 0 && v3 < 0 && v4 < 0 &&
                    v5 < 0 && v6 < 0 && v7 < 0 && v8 < 0)
                  volumeFraction += 1;
                else if (v1 * v2 * v3 * v4 * v5 * v6 * v7 * v8 < 0)
                  volumeFraction += 0.5;
              }
        
        volumeFraction /= NX * NY * NZ;
        //std::cout << "volume fraction is " << volumeFraction << endl;
      }
    
    return (volumeFraction < 0.5);
  }
  
  oRefinementCriteriaOnSign::oRefinementCriteriaOnSign(const oMapping& m, const oLevelSet& ls_)
    : mapping(m), ls(ls_), dim(m.getDim())  {}
  bool oRefinementCriteriaOnSign::operator()(oOctree::const_iterator cell, int level,
                                             const int * ijk,
                                             oOctree::const_iterator children_beg,
                                             oOctree::const_iterator children_end) const
  {
    const bool debug = false;
    std::vector<double> ls_vals(mapping.getTopo().pow_base2[dim]);
    ls.getLevelSetValues(ijk, level, ls_vals);
    
    if (debug) cout << " in coarsening_criteria level is " << level << endl;
    
    double change_sign = *std::max_element(ls_vals.begin(), ls_vals.end()) *
        *std::min_element(ls_vals.begin(), ls_vals.end());
    return (change_sign <= 0.);
  }
  
  oRefinementCriteriaOnSignConstrained::oRefinementCriteriaOnSignConstrained(const oMapping& m, const oLevelSet& ls_, int levelMin_, int targetLevelMax_)
    : mapping(m), ls(ls_), dim(m.getDim()), levelMin(levelMin_), targetLevelMax(targetLevelMax_)  {}
  bool oRefinementCriteriaOnSignConstrained::operator()(oOctree::const_iterator cell, int level,
                                                        const int * ijk,
                                                        oOctree::const_iterator children_beg,
                                                        oOctree::const_iterator children_end) const
  {
    const bool debug = false;
    std::vector<double> ls_vals(mapping.getTopo().pow_base2[dim]);
    ls.getLevelSetValues(ijk, level, ls_vals);
    
    if (debug) cout << " in coarsening_criteria level is " << level << endl;
    
    double change_sign = *std::max_element(ls_vals.begin(), ls_vals.end()) *
        *std::min_element(ls_vals.begin(), ls_vals.end());
    return ((change_sign <= 0. && level < targetLevelMax) || level<levelMin );
  }
  
  
  oRefinementCriteriaPositive::oRefinementCriteriaPositive(const oMapping& m, const oLevelSet& ls_)
    : mapping(m), ls(ls_), dim(m.getDim())  {}
  bool oRefinementCriteriaPositive::operator()(oOctree::const_iterator cell, int level,
                                               const int * ijk,
                                               oOctree::const_iterator children_beg,
                                               oOctree::const_iterator children_end) const
  {
    const bool debug = false;
    std::vector<double> ls_vals(mapping.getTopo().pow_base2[dim]);
    ls.getLevelSetValues(ijk, level, ls_vals);
    
    if (debug) cout << " in coarsening_criteria level is " << level << endl;
    
    double minValue = *std::min_element(ls_vals.begin(), ls_vals.end());
    return (minValue < 0.);
  }

  oRefinementCriteriaOnBand::oRefinementCriteriaOnBand(const oMapping& m,
                                                       const oLevelSet& ls_, const double& d, int levelMin_)
    : mapping(m), ls(ls_), dim(m.getDim()), dist(d), levelMin(levelMin_) {}
  bool oRefinementCriteriaOnBand::operator()(oOctree::const_iterator cell, int level,
                                             const int * ijk,
                                             oOctree::const_iterator children_beg,
                                             oOctree::const_iterator children_end) const
  {
    const bool debug = false;
    std::vector<double> lsv(mapping.getTopo().pow_base2[dim]);
    ls.getLevelSetValues(ijk, level, lsv);

    if (debug) cout << " in coarsening_criteria level is " << level << endl;
    return (min(fabs( *std::max_element(lsv.begin(), lsv.end())),
                fabs( *std::min_element(lsv.begin(), lsv.end()))) <= dist || level < levelMin);

  }


  oRefinementCriteriaOnFunction::oRefinementCriteriaOnFunction(const oMapping& m,
                                                       const oLevelSet& ls_, std::function<bool (std::vector<double> &lsVals,
                                                                                                 std::array<double, 3> &binf,
                                                                                                 std::array<double, 3> &bsup)> func_, int levelMin_)
    : mapping(m), ls(ls_), dim(m.getDim()), levelMin(levelMin_), func(func_) { }




  bool oRefinementCriteriaOnFunction::operator()(oOctree::const_iterator cell, int level,
                                             const int * ijk,
                                             oOctree::const_iterator children_beg,
                                             oOctree::const_iterator children_end) const
  {
    const bool debug = false;
    std::vector<double> lsv(mapping.getTopo().pow_base2[dim]);
    ls.getLevelSetValues(ijk, level, lsv);

    if (debug) cout << " in coarsening_criteria level is " << level << endl;

    std::array<double, 3> binf{0,0,0}, bsup{0,0,0};
    mapping.getBox(ijk, level,&binf[0],&bsup[0]);

    return (func(lsv,binf, bsup) || level < levelMin);
  }





  oRefinementCriteriaOnBandAndPositiveLevelset::oRefinementCriteriaOnBandAndPositiveLevelset(const oMapping& m)
    : mapping(m), dim(m.getDim()){}

  void oRefinementCriteriaOnBandAndPositiveLevelset::addLevelSet(const oLevelSet *ls, double d, int level, int sign){
    ls_list.push_back(ls);
    d_list.push_back(d);
    level_list.push_back(level);
    sign_list.push_back(sign);
  }

  bool oRefinementCriteriaOnBandAndPositiveLevelset::operator()(oOctree::const_iterator cell, int level,
                                                                const int * ijk,
                                                                oOctree::const_iterator children_beg,
                                                                oOctree::const_iterator children_end) const
  {
    bool result=false,res;
    int count=0;
    const bool debug = false;
    std::vector<const oLevelSet* >::const_iterator it = ls_list.begin();
    std::vector<double>::const_iterator it_d=d_list.begin();
    std::vector<int>::const_iterator it_level=level_list.begin();
    std::vector<int>::const_iterator it_sign=sign_list.begin();
    std::vector<const oLevelSet* >::const_iterator it_end = ls_list.end();

    for (;it!=it_end;it++,it_d++,it_level++,it_sign++,count++){
        std::vector<double> lsv(mapping.getTopo().pow_base2[dim]);
        (*it)->getLevelSetValues(ijk, level, lsv);

        if (debug) cout << " in coarsening_criteria level is " << level << endl;

        double maxx =  (*std::max_element(lsv.begin(), lsv.end()));
        double minn =  (*std::min_element(lsv.begin(), lsv.end()));

        if ((*it_sign)>0 && (maxx>=0) && level<(*it_level))// refine if positive levelset... and if necessary according to level...
          res = true;
        else if ((*it_sign)<0 && (minn<=0) && level<(*it_level))// refine if negative levelset... and if necessary according to level...
          res = true;
        else
          res = ((min(fabs(maxx), fabs(minn)) <= (*it_d)) && level<(*it_level));
        if (debug) cout << "distance to levelset " << count << " : " << (*it_d) << " and " << lsv[0] << " " << lsv[1] << " " << lsv[2] << " count = " << count << " res = " << res << endl;

        if (count==0) result=res;
        else result=(result||res);
      }
    if (debug) cout << "Final result " << result << endl;
    return result;
  }






  //It may not be a great idea to refine by curvature everywhere
  //For example in the case of a sphere it will refine heavily near the centre
  //Should be used only elements crossed by the iso-zero
  oRefinementCriteriaOnCurvature::oRefinementCriteriaOnCurvature(const oMapping& m,
                                                                 const oLevelSet& ls_, const double& fact)
    : mapping(m), ls(ls_), dim(m.getDim()), factor(fact) {}
  bool oRefinementCriteriaOnCurvature::operator()(oOctree::const_iterator cell, int level,
                                                  const int * ijk,
                                                  oOctree::const_iterator children_beg,
                                                  oOctree::const_iterator children_end) const
  {


    std::vector<double> lsv(mapping.getTopo().pow_base2[dim]);
    ls.getLevelSetValues(ijk, level, lsv);

    bool refine_by_iso_zero = (*std::max_element(lsv.begin(), lsv.end()) * *std::min_element(lsv.begin(), lsv.end())) <= 0.;

    if(refine_by_iso_zero){

        std::vector<double> lsc(mapping.getTopo().pow_base2[dim]);
        ls.getLevelSetCurvatures(ijk, level, lsc);

        //On veut h < fact * R <=> h/R < fact <=> h * kappa < fact
        const double *step = mapping.getStep(level);
        double h = max(max(step[0], step[1]), step[2]);
        //    cout<<"h="<<h<<endl;
        //    throw;

        bool refine_by_curvature  = (max( fabs( *std::max_element(lsc.begin(), lsc.end())), fabs( *std::min_element(lsc.begin(), lsc.end())))
                                     * h >= factor);

//        cout<<level<<" h="<< h <<" max curv="<<max( fabs( *std::max_element(lsc.begin(), lsc.end())), fabs( *std::min_element(lsc.begin(), lsc.end())))<<
//              " min Rad="<<1/max( fabs( *std::max_element(lsc.begin(), lsc.end())), fabs( *std::min_element(lsc.begin(), lsc.end())));
//        cout<< " crit="<<max( fabs( *std::max_element(lsc.begin(), lsc.end())), fabs( *std::min_element(lsc.begin(), lsc.end()))) * h<<" crit1="<<refine_by_curvature<<endl;

        return refine_by_curvature;
      }


    return false;


  }









} // end namespace


// void oLevelSetOnActive::detect_action(const oOctree::cell_type* cell, int level, const int* ijk)
// {

//   const bool debug = false;
//   bool exist;
//   int ijkn[3];
//   const int* per = octree.getPeriodicity();
//   if (dim == 3)
//     {
//       fill(action_edge,  action_edge+12, 0);
//       fill(action_face,  action_face+6, 0);
//       //loop over the faces
//       int on_what = 2;
//       int count = 0;
//       while(topo.next_neighbor(ijk, level, ijkn, on_what, exist, dim, per, count)) 
// 	{  
// 	  if (!exist || (*octree.cartesian2octree(ijkn, level) > 0))
// 	    {
// 	      //oOctree::cell_type* neighbor = octree.cartesian2octree(ijkn, level);
// 	      //if (*neighbor == 0) 
// 	      //{
// 		  action_face[count] = 1;
// 		  //loop over the edges to set action on the edges
// 		  for (int ed = 0; ed < topo.nb_entities[dim][1]; ++ed)
// 		    action_edge[topo.face_edge_connect[count][ed]] += 1;
// 		  //}
// 	    }
// 	  count++;
// 	}
//       //loop over the edges
//       on_what = 1;
//       count = 0;
//       while(topo.next_neighbor(ijk, level, ijkn, on_what, exist, dim, per, count)) 
// 	{
// 	  if (!exist || (*octree.cartesian2octree(ijkn, level) > 0))
// 	    {
// 	      action_edge[count] += 1;
// 	    }
// 	  count++;
// 	}
//     }
//   if (dim == 2)
//     {
//       fill(action_edge,  action_edge+4, 0);
//       //loop over the edges
//       int on_what = 1;
//       int count = 0;
//       while(topo.next_neighbor(ijk, level, ijkn, on_what, exist, dim, per, count)) 
// 	{
//           if (debug) cout << " exist " << exist << endl;
// 	  if (debug && exist) 
// 	    {
// 	      oOctree::cell_type* neighbor = octree.cartesian2octree(ijkn, level);
// 	      if (*neighbor == 0) cout << " neighbor activity is 0" << endl;
// 	    }
// 	  if (!exist || (*octree.cartesian2octree(ijkn, level) > 0))
// 	    {	
// 	      action_edge[count] = 1;
// 	    }
// 	  count++;
// 	}
//       if (debug) 
// 	{
// 	  cout << "action edge ";
// 	  std::copy(action_edge, action_edge+4, std::ostream_iterator<double>(std::cout, " ")); cout << endl; 
// 	}

//     }
// }
//   oLevelSetOnActive::oLevelSetOnActive(const oOctree& o, const oActiveNodes& active_nodes,  const oLevelSet& lsa) 
//     : oLevelSet(o.getMapping(), o.getLevelMax()), octree(o)
//   {
//     const bool debug = false;    
//     oActiveNodes::const_iterator it  = active_nodes.begin();
//     oActiveNodes::const_iterator ite = active_nodes.end();    
//     for (; it != ite; it++)
//       {
// 	const oKey& key = it->first;
// 	const oKeyInfo& info = it->second;
// 	if (info.isFree()) ls[key] = lsa.getLevelSet(key.ijk);
// 	else
// 	  {
// 	    if (debug) cout << " in oLevelSetOnActive hanging case " << endl;
// 	    double v = 0.;
// 	    oKeyInfo::key_const_iterator itk = info.begin(), itke = info.end(); 
// 	    for ( ; itk != itke; ++itk)
// 	      {
// 		const oKey& key = *itk;
// 		v += lsa.getLevelSet(key.ijk);
// 	      }
// 	    v /= (double) info.size();
// 	    if (debug) cout << " in oLevelSetOnActive v is  " << v << endl;
// 	    ls[key] = v;
// 	  }
//       }
//   }
//   void oLevelSetOnActive::setHangingValues(const oActiveNodes& active_nodes) 
//   {
//     oActiveNodes::const_iterator it  = active_nodes.begin();
//     oActiveNodes::const_iterator ite = active_nodes.end();    
//     for (; it != ite; it++)
//       {
// 	const oKey& key = it->first;
// 	const oKeyInfo& info = it->second;
// 	if (info.isTied())
// 	  {
// 	    oKeyInfo::key_const_iterator itk = info.begin(), itke = info.end(); 
// 	    std::vector<oVal*> ties;
// 	    for ( ; itk != itke; ++itk)
// 	      {
// 		const oKey& keyo = *itk;
// 		ties.push_back(ls.getValPtr(keyo));
// 	      }
// 	    oValManager::iterator it = ls.find(key);
// 	    if (it != ls.end())
// 	      {
// 		delete it->second;
// 		it->second = new oValLC(ties);	    
// 	      }
// 	    else ls.insert(make_pair(key, new oValLC(ties)));
// 	  }
//       }
//   }

