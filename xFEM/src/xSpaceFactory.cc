/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "xSpaceFactory.h"

namespace xfem
{


xSpaceLagrange::lag_degree_t xSpaceFactoryBase::getHierachicalOrder()
{
    switch (order)
    {
        case 0: return xSpaceLagrange::DEGREE_ZERO;
        case 1: return xSpaceLagrange::DEGREE_ONE;
        case 2: return xSpaceLagrange::DEGREE_TWO;
        case 3: return xSpaceLagrange::DEGREE_THREE;
        default :
                {
                    std::cout<<"Hierachical space suport only order 0,1,2 or 3\n Asked for order "<<order<<std::endl;
                    throw -1;
                }
    }
}

void  xSpaceFactoryBase::error1()
{
    std::cout<<"Using getSpace methode with string argument is only available with scalar physics type for now"<<std::endl;
    throw -1;
}

void xSpaceFactoryBase::setPhysString(std::string s1)
{
    const int n = phys_strings.size();
    if (n>0)
    {
        if (n>1)
            phys_strings.resize(1);
    }
    else
        phys_strings.push_back(s1);
}
void xSpaceFactoryBase::forceSetPhysString(std::string &s1)
{
  //const int n = phys_strings.size();
    phys_strings.resize(1);
    phys_strings[0]=s1;
}

void xSpaceFactoryBase::setPhysString(std::string s1,std::string s2)
{
    const int n = phys_strings.size();
    if (n>0)
    {
        if (n>2)
            phys_strings.resize(2);
        else if (n<2)
            phys_strings.push_back(s2);
    }
    else
    {
        phys_strings.push_back(s1);
        phys_strings.push_back(s2);
    }
}

void xSpaceFactoryBase::setPhysString(std::string s1, std::string s2, std::string s3)
{
    const int n = phys_strings.size();
    if (n>0)
    {
        if (n>3)
            phys_strings.resize(3);
        else if (n<2)
        {
            phys_strings.push_back(s2);
            phys_strings.push_back(s3);
        }
        else if (n<3)
            phys_strings.push_back(s3);
    }
    else
    {
        phys_strings.push_back(s1);
        phys_strings.push_back(s2);
        phys_strings.push_back(s3);
    }
}

void  xSpaceFactoryBase::setSpaceProductOrder(const int order_)
{
    order=order_;
    return;
}
void  xSpaceFactoryBase::setSpaceProductPhysStrings(std::vector<std::string> & phys_strings_ )
{
    phys_strings.resize(phys_strings_.size());
    std::copy(phys_strings_.begin(),phys_strings_.end(),phys_strings.begin());
    return;
}

void  xSpaceFactoryBase::setSpaceProductGenerator(xNonLocalInfoGeneratorForKeysAndFcts * generator_)
{
    generator=generator_;
    return;
}

const std::string & xSpaceFactoryBase::getPhysString( unsigned int id) const
{
    return phys_strings[id];
}

// Hierarchical
template <>
xSpaceLagrange xSpaceFactoryBase::createSpaceProduct<xSpaceLagrange>(std::string &phys ,xSpace::TensorialType_t space_tensorial)
{
    xSpaceLagrange hierachical(phys,space_tensorial,getHierachicalOrder());
    return hierachical;
}
// polynomial Lagrange
template <>
xSpacePolynomialLagrange xSpaceFactoryBase::createSpaceProduct<xSpacePolynomialLagrange>(std::string &phys ,xSpace::TensorialType_t space_tensorial)
{
    assert(order>-1);
    xSpacePolynomialLagrange lag(phys,space_tensorial,order);
    return lag;
}
// polynomial Bernstein
template <>
xSpacePolynomialBernstein xSpaceFactoryBase::createSpaceProduct<xSpacePolynomialBernstein>(std::string &phys ,xSpace::TensorialType_t space_tensorial)
{
    assert(order>-1);
    xSpacePolynomialBernstein bern(phys,space_tensorial,order);
    return bern;
}
// polynomial Octree Lagrange
template <>
xSpacePolynomialOctreeLagrange xSpaceFactoryBase::createSpaceProduct<xSpacePolynomialOctreeLagrange>(std::string &phys ,xSpace::TensorialType_t space_tensorial)
{
    assert(order>-1);
    assert(generator);
    xSpacePolynomialOctreeLagrange lag(phys,space_tensorial,order,generator);
    return lag;
}
// polynomial LagrangeQH
template <>
xSpacePolynomialLagrangeQH xSpaceFactoryBase::createSpaceProduct<xSpacePolynomialLagrangeQH>(std::string &phys ,xSpace::TensorialType_t space_tensorial)
{
    assert(order>-1);
    xSpacePolynomialLagrangeQH lag(phys,space_tensorial,order);
    return lag;
}
// polynomial BernsteinQH
template <>
xSpacePolynomialBernsteinQH xSpaceFactoryBase::createSpaceProduct<xSpacePolynomialBernsteinQH>(std::string &phys ,xSpace::TensorialType_t space_tensorial)
{
    assert(order>-1);
    xSpacePolynomialBernsteinQH bern(phys,space_tensorial,order);
    return bern;
}
// polynomial Legendre
template <>
xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH xSpaceFactoryBase::createSpaceProduct<xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH>(std::string &phys ,xSpace::TensorialType_t space_tensorial)
{
    assert(order>-1);
    xSpacePolynomialHierarchicalLegendreTensorProductSpaceQH lag(phys,space_tensorial,order);
    return lag;
}

template <>
xSpacePolynomialHierarchicalLegendreTrunkSpaceQH xSpaceFactoryBase::createSpaceProduct<xSpacePolynomialHierarchicalLegendreTrunkSpaceQH>(std::string &phys ,xSpace::TensorialType_t space_tensorial)
{
    assert(order>-1);
    xSpacePolynomialHierarchicalLegendreTrunkSpaceQH lag(phys,space_tensorial,order);
    return lag;
}

} // end of namespace
