/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE.md & LICENSE.md files for terms
   and conditions.
 */

#include <iostream>

// xfEM
#include "xStateOfValue.h"
#include "xValue.h"

namespace xfem
{
using std::endl;
xValueError::xValueError() : abs2_loc(0.0), abs2_exa_loc(0.), eng_loc(0.0), eng_exa_loc(0.), vol_loc(0.0) {}

double xValueError::vol_tot = 0.;
double xValueError::abs2_tot = 0.;
double xValueError::eng_tot = 0.;
double xValueError::abs2_exa_tot = 0.;
double xValueError::eng_exa_tot = 0.;
bool xValueError::reduced = false;
xValueError::choice_t xValueError::_choice = xValueError::ABS2;

void xValueError::clear()
{
   vol_tot = 0.;
   abs2_tot = 0.;
   eng_tot = 0.;
   abs2_exa_tot = 0.;
   eng_exa_tot = 0.;
   reduced = false;
}

void xValueError::finalizeAgainABS(MPI_Comm world)
{
   assert(reduced);  // you must call this function only after at least one call to finalize
   // you may call this function many time : it accumulate nicely the new summed term removing extras previous
   // summation
   // It only works if previous calls where made with same communicator (or if some specific tricky summation on different
   // communicator are done knowing how this function was implemented)
   int nb_proc_m1;
   MPI_Comm_size(world, &nb_proc_m1);
   --nb_proc_m1;
   double v[] = {abs2_tot, abs2_exa_tot};
   MPI_Allreduce(MPI_IN_PLACE, v, 2, MPI_DOUBLE, MPI_SUM, world);
   abs2_tot = v[0] - nb_proc_m1 * abs2_tot;
   abs2_exa_tot = v[1] - nb_proc_m1 * abs2_exa_tot;
   reduced = true;
}

void xValueError::finalize(MPI_Comm world)
{
   assert(!reduced);  // you must not call this function many time : group all computation before calling this method
   double v[] = {abs2_tot, abs2_exa_tot, eng_tot, eng_exa_tot, vol_tot};
   MPI_Allreduce(MPI_IN_PLACE, v, 5, MPI_DOUBLE, MPI_SUM, world);
   abs2_tot = v[0];
   abs2_exa_tot = v[1];
   eng_tot = v[2];
   eng_exa_tot = v[3];
   vol_tot = v[4];
   reduced = true;
}
double xValueError::total(const std::string &s)
{
   assert(reduced);  // you must call finalize before calling this function
   if (s == "VOL")
      return vol_tot;
   else if (s == "ABS")
      return std::sqrt(abs2_tot);
   else if (s == "ABS_EXA")
      return std::sqrt(abs2_exa_tot);
   else if (s == "ENG")
      return eng_tot;
   else if (s == "ENG_EXA")
      return eng_exa_tot;
   else if (s == "REL")
      return std::sqrt(abs2_tot / eng_tot);
   else if (s == "REL_EXA")
      return std::sqrt(abs2_exa_tot / eng_exa_tot);
   else if (s == "EFF_EXA")
      return std::sqrt(abs2_tot / abs2_exa_tot);
   else
   {
      assert(0), throw;
   };
}

void xValueError::choice(const std::string &s)
{
   if (s == "ABS2")
      _choice = ABS2;
   else if (s == "ABS2_EXA")
      _choice = ABS2_EXA;
   else if (s == "ABS")
      _choice = ABS;
   else if (s == "ABS_EXA")
      _choice = ABS_EXA;
   else if (s == "VOL")
      _choice = VOL;
   else if (s == "REL")
      _choice = REL;
   else if (s == "REL_EXA")
      _choice = REL_EXA;
   else if (s == "DNS")
      _choice = DNS;
   else if (s == "ENG")
      _choice = ENG;
   else if (s == "ENG_EXA")
      _choice = ENG_EXA;
   else if (s == "EFF_EXA")
      _choice = EFF_EXA;
   else
      assert(0);
}

void xValueError::setVal(double in)
{
   switch (_choice)
   {
      case ABS2:
         abs2_loc = in;
         abs2_tot += in;
         break;
      case VOL:
         vol_loc = in;
         vol_tot += in;
         break;
      case ENG:
         eng_loc = in;
         eng_tot += in;
         break;
      case ABS2_EXA:
         abs2_exa_loc = in;
         abs2_exa_tot += in;
         break;
      case ENG_EXA:
         eng_exa_loc = in;
         eng_exa_tot += in;
         break;
      default:
         assert(0);
         break;
   }
}

double xValueError::getVal() const
{
   assert(reduced);  // you must call finalize befor calling this function
   switch (_choice)
   {
      case ABS2:
         return abs2_loc;
         break;
      case ABS2_EXA:
         return abs2_exa_loc;
         break;
      case ABS:
         return std::sqrt(abs2_loc);
         break;
      case ABS_EXA:
         return std::sqrt(abs2_exa_loc);
         break;
      case VOL:
         return vol_loc;
         break;
      case ENG:
         return eng_loc;
         break;
      case ENG_EXA:
         return eng_exa_loc;
         break;
      case REL:
         return std::sqrt(abs2_loc / eng_tot);
         break;
      case REL_EXA:
         return std::sqrt(abs2_exa_loc / eng_exa_tot);
         break;
      case DNS:
         return std::sqrt(abs2_loc * vol_tot / (eng_tot * vol_loc));
         break;
      case EFF_EXA:
         return std::sqrt(abs2_loc / abs2_exa_loc);
         break;
      default:
         assert(0);
         throw;
         break;
   }
}

//----------------------------------------------------------------------------------------------------------------

template <>
xValue<double> *xCloneValue(const xSingleValue<double> *val, const std::map<xValue<double> *, xValue<double> *> &coresp)
{
   xValue<double> *new_val = new xSingleValue<double>;
   const xStateOfValue *state = val->getState();
   if (state) xCloneState(state, new_val);
   new_val->setVal(val->getVal());
   return new_val;
}

}  // namespace xfem
