/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _exlibris_macro__H
#define  _exlibris_macro__H

#ifdef __GNUC__
#define EXLIBRIS_MACRO_WARNUNUSEDTYPE  __attribute__((unused))
#else
#define EXLIBRIS_MACRO_WARNUNUSEDTYPE
#endif

#endif
