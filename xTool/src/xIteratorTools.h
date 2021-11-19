/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/
#ifndef _xIteratorTools_
#define _xIteratorTools_
#include <boost/iterator/filter_iterator.hpp>
namespace xtool
{
//! xRange is a class that describe a range from begin to end ( [begin, end [ )
template <class ITERATOR>
class xRange
{
  public:
   using iterator = ITERATOR;
   xRange(const ITERATOR &_begin, const ITERATOR &_end) : b(_begin), e(_end) {}
   inline const ITERATOR &begin() const { return b; }
   inline const ITERATOR &end() const { return e; }
   inline bool empty() const { return b == e; }
   inline size_t size() const
   {
      using std::distance;
      return distance(b, e);
   }

  private:
   ITERATOR b;
   ITERATOR e;
};
template <class ITERATOR>
xRange<ITERATOR> make_range(const ITERATOR &_begin, const ITERATOR &_end)
{
   return xRange<ITERATOR>(_begin, _end);
}

template <class ITERATOR>
const ITERATOR &begin(const xRange<ITERATOR> &r)
{
   return r.begin();
}
template <class ITERATOR>
const ITERATOR &end(const xRange<ITERATOR> &r)
{
   return r.end();
}

template <class Predicate, class ITERATOR>
boost::filter_iterator<Predicate, ITERATOR> make_filter_iter(const Predicate &p, const ITERATOR &b, const ITERATOR &e)
{
   return boost::filter_iterator<Predicate, ITERATOR>(p, b, e);
}

template <class PRED, class RANGE>
xRange<boost::filter_iterator<PRED, typename RANGE::iterator>> make_filter_range(const PRED &p, const RANGE &r)
{
   auto itb = make_filter_iter(p, r.begin(), r.end());
   auto ite = make_filter_iter(p, r.end(), r.end());
   return make_range(itb, ite);
}

template <class PRED, class ITER>
xRange<boost::filter_iterator<PRED, ITER>> make_filter_range(const PRED &p, const ITER &b, const ITER &e)
{
   auto itb = make_filter_iter(p, b, e);
   auto ite = make_filter_iter(p, e, e);
   return make_range(itb, ite);
}

template <class ITERATOR, class NEWVALUETYPE>
class convertValueIterator;

template <class ITERATOR, class NEWVALUETYPE>
size_t distance(const convertValueIterator<ITERATOR, NEWVALUETYPE> &b, const convertValueIterator<ITERATOR, NEWVALUETYPE> &e)
{
   return distance(b.it, e.it);
};

/// A class of adaptor of forward iterator.
/// it just rewire the * operator so that it return an object of type NEWVALUETYPE
/// instead of an object of type ITERATOR::value_type. To do that,
/// operator * call a NEWVALUETYPE constructor on the original  ITERATOR::value_type object.
template <class ITERATOR, class NEWVALUETYPE>
class convertValueIterator
{
  public:
   typedef convertValueIterator<ITERATOR, NEWVALUETYPE> iter_t;
   typedef NEWVALUETYPE value_type;
   convertValueIterator(const ITERATOR &_it) : it{_it} {}
   NEWVALUETYPE operator*() { return NEWVALUETYPE((*it)); }
   iter_t &operator++()
   {
      ++it;
      return *this;
   }
   bool operator!=(const iter_t &in) { return it != in.it; }
   friend size_t distance<>(const iter_t &b, const iter_t &e);

  private:
   ITERATOR it;
};

/// A adaptor class of forward iterator casting value pointed by iterator ITERATOR
/// by CAST when  dereferencing
template <class ITERATOR, class CAST>
class castValueForwardIterator
{
   static_assert(std::is_same<typename ITERATOR::iterator_category, std::forward_iterator_tag>::value,
                 "ITERATOR type do not correspond to a forward iterator");

  public:
   typedef CAST value_type;
   typedef typename ITERATOR::difference_type difference_type;
   typedef CAST *pointer;
   typedef CAST &reference;
   typedef typename ITERATOR::iterator_category iterator_category;
   typedef castValueForwardIterator<ITERATOR, CAST> iter_t;
   castValueForwardIterator(const ITERATOR &_it) : it{_it} {}
   void operator=(ITERATOR &_it) { it = _it; }
   iter_t &operator++()
   {
      ++it;
      return *this;
   }
   iter_t operator++(int)
   {
      iter_t rit(it);
      ++it;
      return rit;
   }
   value_type operator*() { return static_cast<value_type>(*it); }
   bool operator!=(const iter_t &rhs) { return it != rhs.it; }
   bool operator==(const iter_t &rhs) { return it == rhs.it; }

  private:
   ITERATOR it;
};
/// A adaptor class of bidirectional iterator casting value pointed by iterator ITERATOR
/// by CAST when  dereferencing
template <class ITERATOR, class CAST>
class castValueBidirIterator
{
   static_assert(std::is_same<typename ITERATOR::iterator_category, std::bidirectional_iterator_tag>::value,
                 "ITERATOR type do not correspond to a bidirectional iterator");

  public:
   typedef CAST value_type;
   typedef typename ITERATOR::difference_type difference_type;
   typedef CAST *pointer;
   typedef CAST &reference;
   typedef typename ITERATOR::iterator_category iterator_category;
   typedef castValueBidirIterator<ITERATOR, CAST> iter_t;
   castValueBidirIterator(const ITERATOR &_it) : it{_it} {}
   void operator=(ITERATOR &_it) { it = _it; }
   iter_t &operator++()
   {
      ++it;
      return *this;
   }
   iter_t &operator--()
   {
      ++it;
      return *this;
   }
   iter_t operator++(int)
   {
      iter_t rit(it);
      ++it;
      return rit;
   }
   iter_t operator--(int)
   {
      iter_t rit(it);
      --it;
      return rit;
   }
   value_type operator*() { return static_cast<value_type>(*it); }
   bool operator!=(const iter_t &rhs) { return it != rhs.it; }
   bool operator==(const iter_t &rhs) { return it == rhs.it; }

  protected:
   ITERATOR it;
};

template <typename ITERATOR, typename CAST, typename CAT>
struct castValueIteratorTraits
{
};
template <typename ITERATOR, typename CAST>
struct castValueIteratorTraits<ITERATOR, CAST, std::bidirectional_iterator_tag>
{
   typedef castValueBidirIterator<ITERATOR, CAST> classtype;
};
template <typename ITERATOR, typename CAST>
struct castValueIteratorTraits<ITERATOR, CAST, std::forward_iterator_tag>
{
   typedef castValueForwardIterator<ITERATOR, CAST> classtype;
};

/// A function to create adptor class of iterator casting value pointed by iterator ITERATOR
/// by CAST when  dereferencing. For now only instance of castValueForwardIterator and castValueBidirIterator
/// are created by thys function.
template <typename CAST, typename ITERATOR>
inline auto makeCastValueIterator(const ITERATOR &_it) ->
    typename castValueIteratorTraits<ITERATOR, CAST, typename ITERATOR::iterator_category>::classtype
{
   typedef typename castValueIteratorTraits<ITERATOR, CAST, typename ITERATOR::iterator_category>::classtype classtype;
   return classtype(_it);
}

}  // namespace xtool

#endif
