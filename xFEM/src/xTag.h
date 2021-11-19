/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef XFEM_EXCHANGE_TAG_HH
#define XFEM_EXCHANGE_TAG_HH

// ----------------------------------------------------------------------------
// HEADERS
// ----------------------------------------------------------------------------


namespace xfem
{
// ----------------------------------------------------------------------------
// CLASS Tag
// ----------------------------------------------------------------------------
/*! \ingroup Parallel
    \brief Tag object used to label data exchange between partitions.

    Simple integer encapsulation.

    \todo Other tag design?.
*/
template< typename T >
class Tag
{
  public:
    Tag( int ref ) { tag = 500 + ref; }
    int operator()() const { return tag; }
  private:
    int tag;

};

} // end of namespace
#endif
// == END OF FILE =============================================================
