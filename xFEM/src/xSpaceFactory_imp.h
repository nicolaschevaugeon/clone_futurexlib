/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef __SPACEFACTORYIMP_H
#define __SPACEFACTORYIMP_H


namespace xfem
{

// SCALAR
template <typename SPACE>
class xSpaceFactory<SPACE,xSpaceFactoryBase::SCALAR> : public xSpaceFactoryBase
{
    public:
        xSpaceFactory(){generator=nullptr;order=-999;}
        xSpace::spacePtr  getSpace() override; 
        xSpace::spacePtr  getSpace(std::string s) override; 
    private:
};
template <typename SPACE>
xSpace::spacePtr  xSpaceFactory<SPACE,xSpaceFactoryBase::SCALAR>::getSpace()
{
    setPhysString("scalar");
    return xSpace::spacePtr( new SPACE (createSpaceProduct<SPACE>(phys_strings[0], xSpace::SCALAR)));
}
template <typename SPACE>
xSpace::spacePtr  xSpaceFactory<SPACE,xSpaceFactoryBase::SCALAR>::getSpace(std::string s)
{
    forceSetPhysString(s);
    return xSpace::spacePtr( new SPACE (createSpaceProduct<SPACE>(phys_strings[0], xSpace::SCALAR)));
}
// V1Dx
template <typename SPACE>
class xSpaceFactory<SPACE,xSpaceFactoryBase::V1Dx> : public xSpaceFactoryBase
{
    public:
        xSpaceFactory(){generator=NULL;order=-999;}
        xSpace::spacePtr  getSpace() override; 
        xSpace::spacePtr  getSpace(std::string s) override { error1(); throw;}
    private:
};
template <typename SPACE>
xSpace::spacePtr  xSpaceFactory<SPACE,xSpaceFactoryBase::V1Dx>::getSpace()
{
    setPhysString("DISPLACEMENT_X");
    return xSpace::spacePtr( new SPACE (createSpaceProduct<SPACE>(phys_strings[0], xSpace::VECTOR_X)) );
}
// V1Dy
template <typename SPACE>
class xSpaceFactory<SPACE,xSpaceFactoryBase::V1Dy> : public xSpaceFactoryBase
{
    public:
        xSpaceFactory(){generator=NULL;order=-999;}
        xSpace::spacePtr  getSpace() override; 
        xSpace::spacePtr  getSpace(std::string s) override { error1(); throw;}
    private:
};
template <typename SPACE>
xSpace::spacePtr  xSpaceFactory<SPACE,xSpaceFactoryBase::V1Dy>::getSpace()
{
    setPhysString("DISPLACEMENT_Y");
    return xSpace::spacePtr( new SPACE (createSpaceProduct<SPACE>(phys_strings[0], xSpace::VECTOR_Y)) );
}
// V1Dz
template <typename SPACE>
class xSpaceFactory<SPACE,xSpaceFactoryBase::V1Dz> : public xSpaceFactoryBase
{
    public:
        xSpaceFactory(){generator=NULL;order=-999;}
        xSpace::spacePtr  getSpace() override; 
        xSpace::spacePtr  getSpace(std::string s) override { error1(); throw;}
    private:
};
template <typename SPACE>
xSpace::spacePtr  xSpaceFactory<SPACE,xSpaceFactoryBase::V1Dz>::getSpace()
{
    setPhysString("DISPLACEMENT_Z");
    return xSpace::spacePtr( new SPACE (createSpaceProduct<SPACE>(phys_strings[0], xSpace::VECTOR_Z)) );
}
// V2Dxy
template <typename SPACE>
class xSpaceFactory<SPACE,xSpaceFactoryBase::V2Dxy> : public xSpaceFactoryBase
{
    public:
        xSpaceFactory(){generator=nullptr;order=-999;}
        xSpace::spacePtr  getSpace() override; 
        xSpace::spacePtr  getSpace(std::string s) override { error1(); throw;}
    private:
};
template <typename SPACE>
xSpace::spacePtr  xSpaceFactory<SPACE,xSpaceFactoryBase::V2Dxy>::getSpace()
{
    setPhysString("DISPLACEMENT_X","DISPLACEMENT_Y");
    SPACE  X( createSpaceProduct<SPACE>(phys_strings[0], xSpace::VECTOR_X) );
    SPACE  Y( createSpaceProduct<SPACE>(phys_strings[1], xSpace::VECTOR_Y) );
    return xSpace::spacePtr(new xSpaceComposite(X,Y));
}
// V2Dxz
template <typename SPACE>
class xSpaceFactory<SPACE,xSpaceFactoryBase::V2Dxz> : public xSpaceFactoryBase
{
    public:
        xSpaceFactory(){generator=NULL;order=-999;}
        xSpace::spacePtr  getSpace() override; 
        xSpace::spacePtr  getSpace(std::string s) override { error1(); throw;}
    private:
};
template <typename SPACE>
xSpace::spacePtr  xSpaceFactory<SPACE,xSpaceFactoryBase::V2Dxz>::getSpace()
{
    setPhysString("DISPLACEMENT_X","DISPLACEMENT_Z");
    SPACE  X( createSpaceProduct<SPACE>(phys_strings[0], xSpace::VECTOR_X) );
    SPACE  Z( createSpaceProduct<SPACE>(phys_strings[1], xSpace::VECTOR_Z) );
    return xSpace::spacePtr(new xSpaceComposite(X,Z));
}
// V2Dyz
template <typename SPACE>
class xSpaceFactory<SPACE,xSpaceFactoryBase::V2Dyz> : public xSpaceFactoryBase
{
    public:
        xSpaceFactory(){generator=NULL;order=-999;}
        xSpace::spacePtr  getSpace() override; 
        xSpace::spacePtr  getSpace(std::string s) override { error1(); throw;}
    private:
};
template <typename SPACE>
xSpace::spacePtr  xSpaceFactory<SPACE,xSpaceFactoryBase::V2Dyz>::getSpace()
{
    setPhysString("DISPLACEMENT_Y","DISPLACEMENT_Z");
    SPACE  Y( createSpaceProduct<SPACE>(phys_strings[0], xSpace::VECTOR_Y) );
    SPACE  Z( createSpaceProduct<SPACE>(phys_strings[1], xSpace::VECTOR_Z) );
    return xSpace::spacePtr(new xSpaceComposite(Y,Z));
}
// V3D
template <typename SPACE>
class xSpaceFactory<SPACE,xSpaceFactoryBase::V3D> : public xSpaceFactoryBase
{
    public:
        xSpaceFactory(){generator=nullptr;order=-999;}
        xSpace::spacePtr  getSpace() override; 
        xSpace::spacePtr  getSpace(std::string s) override { error1(); throw;}
    private:
};
template <typename SPACE>
xSpace::spacePtr  xSpaceFactory<SPACE,xSpaceFactoryBase::V3D>::getSpace()
{
    setPhysString("DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z");
    SPACE  X( createSpaceProduct<SPACE>(phys_strings[0], xSpace::VECTOR_X) );
    SPACE  Y( createSpaceProduct<SPACE>(phys_strings[1], xSpace::VECTOR_Y) );
    SPACE  Z( createSpaceProduct<SPACE>(phys_strings[2], xSpace::VECTOR_Z) );
    return xSpace::spacePtr(new xSpaceComposite(X,Y,Z));
}

} // end of namespace

#endif














