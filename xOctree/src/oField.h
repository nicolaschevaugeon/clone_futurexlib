/*
    octree is a subproject of  xfem : C++ Finite Element Library
    developed under the GNU Lesser General Public License
    See the NOTICE, CONTRIBUTORS & LICENSE files for conditions.
*/

// -*- C++ -*-

#ifndef _OFIELD_H__
#define _OFIELD_H__


namespace xoctree {




class oField
{


public:
  typedef std::function<double (const int*)> f_ijk_t;
  struct oZeroVal
  {
    double operator()(const int*) const { return 0; }
  };

  oField(oKeyManager& key_manager_, f_ijk_t f = oZeroVal()) : key_manager(key_manager_)
  {
    field_name = key_manager.registerField(f);
  }

  double getVal(const oKey* key) const
  {
    oKey::setField(field_name);
    return key->getVal();
  }
  void setVal(oKey* key, const double& v) const
  {
    oKey::setField(field_name);
    return key->setVal(v);
  }
  const oKeyManager& getKeyManager() const {return key_manager;}


private:
  oKeyManager& key_manager;
  int field_name;
};


/// Cell data class
/// Does not work for aniso octree for the moment !!!!
class oCellFieldOnFine{

public:
    using oCellFieldContainer = std::vector<double>;
    using ijkToLinear_t = std::function<int (const int*)>;

    oCellFieldOnFine(oOctree &octree_) : octree(octree_),
        powp(octree_.getTopo().pow_base2[octree_.getLevelMax()]),
        ijkToLinear([this](const int *ijk){return ijk[0] + ijk[1] * powp + ijk[2] * powp * powp;})
    {

        if(octree.getDim() == 2) container = new oCellFieldContainer(powp * powp);
        else container = new oCellFieldContainer(powp * powp * powp);

    }


    oCellFieldOnFine(oOctree &octree_, oCellFieldContainer &containerIn) : octree(octree_),
        powp(octree_.getTopo().pow_base2[octree_.getLevelMax()]),
        ijkToLinear([this](const int *ijk){return ijk[0] + ijk[1] * powp + ijk[2] * powp * powp;}),
        container(&containerIn)
    {

        assert(container->size() == static_cast<size_t>((octree.getDim() == 2)?powp * powp:powp * powp * powp));
        /*
        if(octree.getDim() == 2) assert(container->size() == powp * powp);
        else assert( container->size() == powp * powp * powp);
        */

    }


    void setVal(const int *ijk, const double v){
        (*container)[ijkToLinear(ijk)] = v;
        return;
    }

    double getVal(const int *ijk) const {
        return (*container)[ijkToLinear(ijk)];
    }

    int size() const {
        return container->size();
    }


private:
    oOctree &octree;
    const int powp;
    const ijkToLinear_t ijkToLinear;
    oCellFieldContainer *container{nullptr};
};



} //end namespace 

#endif





