/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _xDATAEXCHANGERTRAITS_H
#define _xDATAEXCHANGERTRAITS_H

namespace xtool
{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Data key association traits ////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Note : the keys_association traits are not used so far in xDataExchanger.h, but for some type def in xKeyInformationContainer
struct send_and_recv_keys_association_trait {};
struct send_or_recv_keys_association_trait {};
struct no_keys_association_trait {};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Data style traits //////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Your info manager need one of the folling traits :
///   typedef xxx_data_style_trait data_style_trait
struct homogeneous_data_style_trait {};
struct nonhomogeneous_data_style_trait {};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// communication pattern traits /////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Your info manager need one of the folling traits :
///   typedef xxx_communication_trait communication_trait
struct send_and_recv_keys_communication_trait {};
struct send_only_keys_communication_trait {};
struct recv_only_keys_communication_trait {};

} // end of namespace

#endif
