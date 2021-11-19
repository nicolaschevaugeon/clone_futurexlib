/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef  _xFactory_H
#define  _xFactory_H

namespace xfem  {

////////////////////////////////////////////////////////////////////////////////
// class template DefaultFactoryError
// Manages the "Unknown Type" error in an object factory
////////////////////////////////////////////////////////////////////////////////

    template <typename IdentifierType, class AbstractProduct>
    struct xDefaultFactoryError
    {
        struct xException : public std::exception
        {
            const char* what() const throw() override { return "Unknown Type in xFactory.h"; }
        };
        
        static AbstractProduct* onUnknownType(IdentifierType)
        {
            throw xException();
        }
    };

////////////////////////////////////////////////////////////////////////////////
// class template Factory
// Implements a generic object factory
////////////////////////////////////////////////////////////////////////////////

    template
    <
        class AbstractProduct, 
        typename IdentifierType,
        typename ProductCreator = AbstractProduct* (*)(),
        template<typename, class>
            class FactoryErrorPolicy = xDefaultFactoryError
    >
    class xFactory 
        : public FactoryErrorPolicy<IdentifierType, AbstractProduct>
    {
    public:
        bool registerCreator(const IdentifierType& id, ProductCreator creator)
        {
            return associations_.insert(
                typename IdToProductMap::value_type(id, creator)).second;
        }
        
        bool unRegisterCreator(const IdentifierType& id)
        {
            return associations_.erase(id) == 1;
        }
        
        AbstractProduct* createObject(const IdentifierType& id)
        {
            typename IdToProductMap::iterator i = associations_.find(id);
            if (i != associations_.end())
            {
                return (i->second)();
            }
            return this->onUnknownType(id);
        }
        
    private:
        typedef std::map<IdentifierType, ProductCreator> IdToProductMap;
        IdToProductMap associations_;
    };


} // namespace xfem

#endif
