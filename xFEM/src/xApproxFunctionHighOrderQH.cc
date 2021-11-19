/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#include "xApproxFunctionHighOrderQH.h"

#include <cmath>

#include "xApproxFunctionHighOrder.h"
#include "xGeomElem.h"

using AOMD::mEntity;
using std::cout;
using std::endl;

namespace xfem
{
int tensorSpaceDefinition::nbShapeFunctionEdgeTotalQH(const int& order) { return order + 1; }

int tensorSpaceDefinition::nbShapeFunctionEdgeQH(const int& order) { return std::max(order - 1, 0); }

int tensorSpaceDefinition::nbShapeFunctionQuadTotal(const int& order) { return (order + 1) * (order + 1); }

int tensorSpaceDefinition::nbShapeFunctionQuad(const int& order) { return (order > 1) ? (order - 1) * (order - 1) : 0; }

int tensorSpaceDefinition::nbShapeFunctionHexTotal(const int& order) { return (order + 1) * (order + 1) * (order + 1); }

int tensorSpaceDefinition::nbShapeFunctionHex(const int& order)
{
   return (order > 1) ? (order - 1) * (order - 1) * (order - 1) : 0;
}

bool tensorSpaceDefinition::acceptShapeFunc(const int& order, const int i, const int j, const int k) { return true; }

// const int xAOMDSwapTable::Hfv[6][4] = {{0,3,2,1},
//                                       {0,1,5,4},
//                                       {1,2,6,5},
//                                       {2,3,7,6},
//                                       {3,0,4,7},
//                                       {4,5,6,7}
//                                      };

// const int xAOMDSwapTable::Hev[12][2] = {{0,1},{1,2},{2,3},{3,0},
//                                        {0,4},{1,5},{2,6},{3,7},
//                                        {4,5},{5,6},{6,7},{7,4}
//                                       };

// const int xAOMDSwapTable::Fev[4][2] = {{0,1},{1,2},{2,3},{3,0}};

const int xSzaboSwapTable::Hfv[6][4] = {{0, 1, 2, 3}, {0, 1, 5, 4}, {1, 2, 6, 5}, {3, 2, 6, 7}, {0, 3, 7, 4}, {4, 5, 6, 7}};

const int xSzaboSwapTable::Hev[12][2] = {{0, 1}, {1, 2}, {3, 2}, {0, 3}, {0, 4}, {1, 5},
                                         {2, 6}, {3, 7}, {4, 5}, {5, 6}, {7, 6}, {4, 7}};

const int xSzaboSwapTable::Fev[4][2] = {{0, 1}, {1, 2}, {3, 2}, {0, 3}};

int trunkSpaceDefinition::nbShapeFunctionEdgeTotalQH(const int& order) { return order + 1; }

int trunkSpaceDefinition::nbShapeFunctionEdgeQH(const int& order) { return std::max(order - 1, 0); }

int trunkSpaceDefinition::nbShapeFunctionQuadTotal(const int& order)
{
   return (order > 3) ? 4 + 4 * (order - 1) + (order - 2) * (order - 3) / 2 : 4 + 4 * (order - 1);
}

int trunkSpaceDefinition::nbShapeFunctionQuad(const int& order) { return (order > 3) ? (order - 2) * (order - 3) / 2 : 0; }

int trunkSpaceDefinition::nbShapeFunctionHexTotal(const int& order)
{
   if (order < 4)
      return 8 + 12 * (order - 1);
   else if (order < 6)
      return 8 + 12 * (order - 1) + 3 * (order - 2) * (order - 3);
   else
      return 8 + 12 * (order - 1) + 3 * (order - 2) * (order - 3) + (order - 3) * (order - 4) * (order - 5) / 6;
}

int trunkSpaceDefinition::nbShapeFunctionHex(const int& order)
{
   return (order > 5) ? (order - 3) * (order - 4) * (order - 5) / 6 : 0;
}

bool trunkSpaceDefinition::acceptShapeFunc(const int& order, const int i, const int j, const int k)
{
   return (i + j + k <= order);
}

}  // end namespace xfem

//------------------------- FF Hierarchiques Legendre VERSION DEPRECIEE ------------------------
using xtensor::xPoint;

// Max number of tabulated Lagrange functions
#define XAPPROXFUNCTIONLEGENDREQHMAXORDER 19

namespace xfem
{
///////////////////////////////////////////////////////////////////////////////////////////////
// xApproxFunction Hierarchical Legendre Implementation     QUADS AND HEXAHEDRA              //
///////////////////////////////////////////////////////////////////////////////////////////////

namespace deprecatedLegendre
{
std::function<double(double)> legPolys[16] = {
    [](double u) -> double { return 1.; },
    [](double u) -> double { return u; },
    [](double u) -> double { return 1.5 * u * u - 0.5; },
    [](double u) -> double { return 2.5 * u * u * u - 1.5 * u; },
    [](double u) -> double { return 4.375 * u * u * u * u - 3.75 * u * u + 0.375; },
    [](double u) -> double { return 7.875 * u * u * u * u * u - 8.75 * u * u * u + 1.875 * u; },
    [](double u) -> double { return 14.4375 * u * u * u * u * u * u - 19.6875 * u * u * u * u + 6.5625 * u * u - 0.3125; },
    [](double u) -> double {
       return 26.8125 * u * u * u * u * u * u * u - 43.3125 * u * u * u * u * u + 19.6875 * u * u * u - 2.1875 * u;
    },
    [](double u) -> double {
       return 50.2734375 * u * u * u * u * u * u * u * u - 93.84375 * u * u * u * u * u * u + 54.140625 * u * u * u * u -
              9.84375 * u * u + 0.2734375;
    },
    [](double u) -> double {
       return 94.9609375 * u * u * u * u * u * u * u * u * u - 201.09375 * u * u * u * u * u * u * u +
              140.765625 * u * u * u * u * u - 36.09375 * u * u * u + 2.4609375 * u;
    },
    [](double u) -> double {
       return 180.42578125 * u * u * u * u * u * u * u * u * u * u - 427.32421875 * u * u * u * u * u * u * u * u +
              351.9140625 * u * u * u * u * u * u - 117.3046875 * u * u * u * u + 13.53515625 * u * u - 0.24609375;
    },
    [](double u) -> double {
       return 344.44921875 * u * u * u * u * u * u * u * u * u * u * u - 902.12890625 * u * u * u * u * u * u * u * u * u +
              854.6484375 * u * u * u * u * u * u * u - 351.9140625 * u * u * u * u * u + 58.65234375 * u * u * u -
              2.70703125 * u;
    },
    [](double u) -> double {
       return 660.1943359375 * u * u * u * u * u * u * u * u * u * u * u * u -
              1894.470703125 * u * u * u * u * u * u * u * u * u * u + 2029.7900390625 * u * u * u * u * u * u * u * u -
              997.08984375 * u * u * u * u * u * u + 219.9462890625 * u * u * u * u - 17.595703125 * u * u + 0.2255859375;
    },
    [](double u) -> double {
       return 1269.6044921875 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              3961.166015625 * u * u * u * u * u * u * u * u * u * u * u + 4736.1767578125 * u * u * u * u * u * u * u * u * u -
              2706.38671875 * u * u * u * u * u * u * u + 747.8173828125 * u * u * u * u * u - 87.978515625 * u * u * u +
              2.9326171875 * u;
    },
    [](double u) -> double {
       return 2448.52294921875 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              8252.42919921875 * u * u * u * u * u * u * u * u * u * u * u * u +
              10893.20654296875 * u * u * u * u * u * u * u * u * u * u - 7104.26513671875 * u * u * u * u * u * u * u * u +
              2368.08837890625 * u * u * u * u * u * u - 373.90869140625 * u * u * u * u + 21.99462890625 * u * u - 0.20947265625;
    },
    [](double u) -> double {
       return 4733.81103515625 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              17139.66064453125 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              24757.28759765625 * u * u * u * u * u * u * u * u * u * u * u -
              18155.34423828125 * u * u * u * u * u * u * u * u * u + 7104.26513671875 * u * u * u * u * u * u * u -
              1420.85302734375 * u * u * u * u * u + 124.63623046875 * u * u * u - 3.14208984375 * u;
    }};

std::function<double(double)> dLegPolys[16] = {
    [](double u) -> double { return 0.; },
    [](double u) -> double { return 1.0000000000000000000; },
    [](double u) -> double { return 3.0 * u; },
    [](double u) -> double { return 7.5 * u * u - 1.5; },
    [](double u) -> double { return 2.5 * u * (7.0 * u * u - 3.0); },
    [](double u) -> double { return 39.375 * u * u * u * u - 26.25 * u * u + 1.875; },
    [](double u) -> double { return 2.625 * u * (33.0 * u * u * u * u - 30.0 * u * u + 5.0); },
    [](double u) -> double { return 187.6875 * u * u * u * u * u * u - 216.5625 * u * u * u * u + 59.0625 * u * u - 2.1875; },
    [](double u) -> double {
       return 0.5625 * u * (715.0 * u * u * u * u * u * u - 1001.0 * u * u * u * u + 385.0 * u * u - 35.0);
    },
    [](double u) -> double {
       return 854.6484375 * u * u * u * u * u * u * u * u - 1407.65625 * u * u * u * u * u * u + 703.828125 * u * u * u * u -
              108.28125 * u * u + 2.4609375;
    },
    [](double u) -> double {
       return 0.4296875 * u *
              (4199.0 * u * u * u * u * u * u * u * u - 7956.0 * u * u * u * u * u * u + 4914.0 * u * u * u * u - 1092.0 * u * u +
               63.0);
    },
    [](double u) -> double {
       return 3788.94140625 * u * u * u * u * u * u * u * u * u * u - 8119.16015625 * u * u * u * u * u * u * u * u +
              5982.5390625 * u * u * u * u * u * u - 1759.5703125 * u * u * u * u + 175.95703125 * u * u - 2.70703125;
    },
    [](double u) -> double {
       return 0.15234375 * u *
              (52003.0 * u * u * u * u * u * u * u * u * u * u - 124355.0 * u * u * u * u * u * u * u * u +
               106590.0 * u * u * u * u * u * u - 39270.0 * u * u * u * u + 5775.0 * u * u - 231.0);
    },
    [](double u) -> double {
       return 16504.8583984375 * u * u * u * u * u * u * u * u * u * u * u * u -
              43572.826171875 * u * u * u * u * u * u * u * u * u * u + 42625.5908203125 * u * u * u * u * u * u * u * u -
              18944.70703125 * u * u * u * u * u * u + 3739.0869140625 * u * u * u * u - 263.935546875 * u * u + 2.9326171875;
    },
    [](double u) -> double {
       return 0.1025390625 * u *
              (334305.0 * u * u * u * u * u * u * u * u * u * u * u * u - 965770.0 * u * u * u * u * u * u * u * u * u * u +
               1062347.0 * u * u * u * u * u * u * u * u - 554268.0 * u * u * u * u * u * u + 138567.0 * u * u * u * u -
               14586.0 * u * u + 429.0);
    },
    [](double u) -> double {
       return 71007.16552734375 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              222815.58837890625 * u * u * u * u * u * u * u * u * u * u * u * u +
              272330.16357421875 * u * u * u * u * u * u * u * u * u * u - 163398.09814453125 * u * u * u * u * u * u * u * u +
              49729.85595703125 * u * u * u * u * u * u - 7104.26513671875 * u * u * u * u + 373.90869140625 * u * u -
              3.14208984375;
    }};

// Version 1 call

// for k in range(2,16):
//     print(N((legendre(k,u) - legendre(k-2,u))/sqrt(4.*k-2.),80))
std::function<double(double)> phiFuncs[20] = {
    [](double u) -> double { return 0.; },
    [](double u) -> double { return 0.; },
    [](double u) -> double { return 0.612372435695795 * u * u - 0.61237243569579458135621052861097268760204315185546875; },
    [](double u) -> double { return 0.790569415042095 * u * u * u - 0.790569415042095 * u; },
    [](double u) -> double {
       return 1.16926793336686 * u * u * u * u - 1.40312152004023 * u * u +
              0.233853586673371360848960875955526717007160186767578125;
    },
    [](double u) -> double {
       return 1.85615530061469 * u * u * u * u * u - 2.65165042944955 * u * u * u + 0.795495128834866 * u;
    },
    [](double u) -> double {
       return 3.07808534238413 * u * u * u * u * u * u - 5.13014223730688 * u * u * u * u + 2.19863238741723 * u * u -
              0.146575492494482151339951769841718487441539764404296875;
    },
    [](double u) -> double {
       return 5.25836387339256 * u * u * u * u * u * u * u - 10.0386946673858 * u * u * u * u * u + 5.57705259299211 * u * u * u -
              0.796721798998873 * u;
    },
    [](double u) -> double {
       return 9.17863192069204 * u * u * u * u * u * u * u * u - 19.7693610599521 * u * u * u * u * u * u +
              13.4791098136037 * u * u * u * u - 2.99535773635638 * u * u +
              0.1069770620127277471755888882398721762001514434814453125;
    },
    [](double u) -> double {
       return 16.2856664250562 * u * u * u * u * u * u * u * u * u - 39.0855994201349 * u * u * u * u * u * u * u +
              31.5691379931859 * u * u * u * u * u - 9.56640545248057 * u * u * u + 0.797200454373381 * u;
    },
    [](double u) -> double {
       return 29.2689266430031 * u * u * u * u * u * u * u * u * u * u - 77.4765705255964 * u * u * u * u * u * u * u * u +
              72.31146582389 * u * u * u * u * u * u - 27.8121022399577 * u * u * u * u + 3.79255939635787 * u * u -
              0.08427909769684148455493044593822560273110866546630859375;
    },
    [](double u) -> double {
       return 53.1496683449504 * u * u * u * u * u * u * u * u * u * u * u -
              153.854303103804 * u * u * u * u * u * u * u * u * u + 162.904556227557 * u * u * u * u * u * u * u -
              76.0221262395266 * u * u * u * u * u + 14.6196396614474 * u * u * u - 0.797434890624405 * u;
    },
    [](double u) -> double {
       return 97.3403443330083 * u * u * u * u * u * u * u * u * u * u * u * u -
              305.926796475169 * u * u * u * u * u * u * u * u * u * u + 362.281732667963 * u * u * u * u * u * u * u * u -
              198.899774798097 * u * u * u * u * u * u + 49.7249436995244 * u * u * u * u - 4.58999480303302 * u * u +
              0.06954537580353058190407722349846153520047664642333984375;
    },
    [](double u) -> double {
       return 179.549189170137 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              608.905945881334 * u * u * u * u * u * u * u * u * u * u * u +
              797.376833892223 * u * u * u * u * u * u * u * u * u - 503.606421405614 * u * u * u * u * u * u * u +
              155.52551249291 * u * u * u * u * u - 20.7367349990547 * u * u * u + 0.797566730732873 * u;
    },
    [](double u) -> double {
       return 333.201769393364 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              1212.85444059184 * u * u * u * u * u * u * u * u * u * u * u * u +
              1740.18245824047 * u * u * u * u * u * u * u * u * u * u - 1242.98747017177 * u * u * u * u * u * u * u * u +
              457.942752168545 * u * u * u * u * u * u - 80.8134268532726 * u * u * u * u + 5.38756179021818 * u * u -
              0.0592039757166832603108019839055486954748630523681640625;
    },
    [](double u) -> double {
       return 621.579840858358 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              2417.25493667139 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              3770.91770120737 * u * u * u * u * u * u * u * u * u * u * u -
              3005.80396473051 * u * u * u * u * u * u * u * u * u + 1288.20169917022 * u * u * u * u * u * u * u -
              284.760375606049 * u * u * u * u * u + 27.9176838829459 * u * u * u - 0.797648110941313 * u;
    },
    [](double u) -> double {
       return 1164.81454265326 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              4819.92224546176 * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              8122.46156179668 * u * u * u * u * u * u * u * u * u * u * u * u -
              7147.76617438108 * u * u * u * u * u * u * u * u * u * u + 3496.18997659944 * u * u * u * u * u * u * u * u -
              932.317327093184 * u * u * u * u * u * u + 122.673332512261 * u * u * u * u - 6.18521004263501;
    },
    [](double u) -> double {
       return 2191.52121718174 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              9614.41566247472 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              17405.4076648249 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              16760.7629364981 * u * u * u * u * u * u * u * u * u * u * u +
              9218.41961507394 * u * u * u * u * u * u * u * u * u - 2885.76614037097 * u * u * u * u * u * u * u +
              480.961023395162 * u * u * u * u * u - 36.1624829620423 * u * u * u + 0.79770183004505 * u;
    },
    [](double u) -> double {
       return 4137.74923101061 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              19184.1100710492 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              37130.5356213856 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              38837.6866844378 * u * u * u * u * u * u * u * u * u * u * u * u +
              23734.141862712 * u * u * u * u * u * u * u * u * u * u - 8544.29107057631 * u * u * u * u * u * u * u * u +
              1733.6242751894 * u * u * u * u * u * u - 176.900436243816 * u * u * u * u + 6.98291195699274 * u * u -
              0.0456399474313250730350688399994396604597568511962890625;
    },
    [](double u) -> double {
       return 7836.92065721015 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              38288.955210941 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              78898.453161939 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              89078.8987312215 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              59897.8801813386 * u * u * u * u * u * u * u * u * u * u * u -
              24402.8400738787 * u * u * u * u * u * u * u * u * u + 5856.68161773088 * u * u * u * u * u * u * u -
              763.914993617072 * u * u * u * u * u + 45.4711305724447 * u * u * u - 0.797739132849908 * u;
    }};

//  for k in range(2,16):
//     print(N(diff(legendre(k,u) - legendre(k-2,u),u)/sqrt(4.*k-2.),80))
std::function<double(double)> dPhiFuncs[20] = {
    [](double u) -> double { return 0.; },
    [](double u) -> double { return 0.; },
    [](double u) -> double { return 1.2247448713915891627124210572219453752040863037109375 * u; },
    [](double u) -> double { return 2.37170824512628 * u * u - 0.79056941504209488069676581289968453347682952880859375; },
    [](double u) -> double { return 4.67707173346743 * u * u * u - 2.80624304008046 * u; },
    [](double u) -> double {
       return 9.28077650307344 * u * u * u * u - 7.95495128834866 * u * u +
              0.79549512883486606096283821898396126925945281982421875;
    },
    [](double u) -> double { return 18.4685120543048 * u * u * u * u * u - 20.5205689492275 * u * u * u + 4.39726477483446 * u; },
    [](double u) -> double {
       return 36.8085471137479 * u * u * u * u * u * u - 50.193473336929 * u * u * u * u + 16.7311577789763 * u * u -
              0.796721798998872632324719234020449221134185791015625;
    },
    [](double u) -> double {
       return 73.4290553655363 * u * u * u * u * u * u * u - 118.616166359713 * u * u * u * u * u + 53.9164392544148 * u * u * u -
              5.99071547271275 * u;
    },
    [](double u) -> double {
       return 146.570997825506 * u * u * u * u * u * u * u * u - 273.599195940944 * u * u * u * u * u * u +
              157.845689965929 * u * u * u * u - 28.6992163574417 * u * u + 0.79720045437338082905398550792597234249114990234375;
    },
    [](double u) -> double {
       return 292.689266430031 * u * u * u * u * u * u * u * u * u - 619.812564204771 * u * u * u * u * u * u * u +
              433.86879494334 * u * u * u * u * u - 111.248408959831 * u * u * u + 7.58511879271573 * u;
    },
    [](double u) -> double {
       return 584.646351794454 * u * u * u * u * u * u * u * u * u * u - 1384.68872793423 * u * u * u * u * u * u * u * u +
              1140.3318935929 * u * u * u * u * u * u - 380.110631197633 * u * u * u * u + 43.8589189843423 * u * u -
              0.7974348906244046464308894428540952503681182861328125;
    },
    [](double u) -> double {
       return 1168.0841319961 * u * u * u * u * u * u * u * u * u * u * u - 3059.26796475169 * u * u * u * u * u * u * u * u * u +
              2898.25386134371 * u * u * u * u * u * u * u - 1193.39864878858 * u * u * u * u * u + 198.899774798097 * u * u * u -
              9.17998960606604 * u;
    },
    [](double u) -> double {
       return 2334.13945921178 * u * u * u * u * u * u * u * u * u * u * u * u -
              6697.96540469467 * u * u * u * u * u * u * u * u * u * u + 7176.39150503 * u * u * u * u * u * u * u * u -
              3525.2449498393 * u * u * u * u * u * u + 777.627562464552 * u * u * u * u - 62.2102049971641 * u * u +
              0.79756673073287343012083283610991202294826507568359375;
    },
    [](double u) -> double {
       return 4664.82477150709 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              14554.2532871021 * u * u * u * u * u * u * u * u * u * u * u +
              17401.8245824047 * u * u * u * u * u * u * u * u * u - 9943.89976137412 * u * u * u * u * u * u * u +
              2747.65651301127 * u * u * u * u * u - 323.253707413091 * u * u * u + 10.7751235804364 * u;
    },
    [](double u) -> double {
       return 9323.69761287537 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              31424.3141767281 * u * u * u * u * u * u * u * u * u * u * u * u +
              41480.0947132811 * u * u * u * u * u * u * u * u * u * u - 27052.2356825746 * u * u * u * u * u * u * u * u +
              9017.41189419154 * u * u * u * u * u * u - 1423.80187803024 * u * u * u * u + 83.7530516488378 * u * u -
              0.79764811094131260471584710103343240916728973388671875;
    },
    [](double u) -> double {
       return 18637.0326824522 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              67478.9114364647 * u * u * u * u * u * u * u * u * u * u * u * u * u +
              97469.5387415601 * u * u * u * u * u * u * u * u * u * u * u -
              71477.6617438108 * u * u * u * u * u * u * u * u * u + 27969.5198127955 * u * u * u * u * u * u * u -
              5593.9039625591 * u * u * u * u * u + 490.693330049044 * u * u * u - 12.37042008527 * u;
    },
    [](double u) -> double {
       return 37255.8606920895 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              144216.234937121 * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              226270.299642724 * u * u * u * u * u * u * u * u * u * u * u * u -
              184368.392301479 * u * u * u * u * u * u * u * u * u * u + 82965.7765356655 * u * u * u * u * u * u * u * u -
              20200.3629825968 * u * u * u * u * u * u + 2404.80511697581 * u * u * u * u - 108.487448886127 * u * u +
              0.797701830045050019890595649485476315021514892578125;
    },
    [](double u) -> double {
       return 74479.4861581911 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              306945.761136787 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              519827.498699398 * u * u * u * u * u * u * u * u * u * u * u * u * u -
              466052.240213253 * u * u * u * u * u * u * u * u * u * u * u + 237341.41862712 * u * u * u * u * u * u * u * u * u -
              68354.3285646105 * u * u * u * u * u * u * u + 10401.7456511364 * u * u * u * u * u - 707.601744975264 * u * u * u +
              13.9658239139855 * u;
    },
    [](double u) -> double {
       return 148901.492486993 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              650912.238585997 * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u * u +
              1183476.79742909 * u * u * u * u * u * u * u * u * u * u * u * u * u * u -
              1158025.68350588 * u * u * u * u * u * u * u * u * u * u * u * u +
              658876.681994724 * u * u * u * u * u * u * u * u * u * u - 219625.560664908 * u * u * u * u * u * u * u * u +
              40996.7713241162 * u * u * u * u * u * u - 3819.57496808536 * u * u * u * u + 136.413391717334 * u * u -
              0.79773913284990793925999241764657199382781982421875;
    }};

inline double phiFunc(int k, double u) { return phiFuncs[k](u); }

inline double dPhiFunc(int k, double u) { return dPhiFuncs[k](u); }

}  // namespace deprecatedLegendre

xApproxFunctionScalarEdgeHierarchicalLegendreQH_deprecated::xApproxFunctionScalarEdgeHierarchicalLegendreQH_deprecated(
    int _order, int _edge, int _ionedge, int _swap)
    : order(_order), ionedge(_ionedge), edge(_edge), swap(_swap)
{
   if (order > XAPPROXFUNCTIONLEGENDREQHMAXORDER)
   {
      cout << "Tabulated functions up to order XAPPROXFUNCTIONLEGENDREQHMAXORDER !\n";
      throw;
   }
}

void xApproxFunctionScalarEdgeHierarchicalLegendreQH_deprecated::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                                        double& res) const
{
   xtensor::xPoint uvw = geo_appro->getUVW();

   // Principe (different des FF lagrange et bernstein):
   // REM : k commence a 2 pour les deprecatedLegendre::phiFuncs, d'ou le terme "ionedge + 2"

   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);

   // On corrige le cas de l'edge (sinon, on a v et w nuls et donc la fonction de forme aussi)
   if (geo_appro->getEntity()->getType() == mEntity::EDGE)
   {
      v = -1.;
      w = -1.;
   }

   switch (edge)
   {
      case 0:  // 2D et 3D ici...
         if (swap) u *= -1;
         res = 0.5 * deprecatedLegendre::phiFunc(ionedge + 2, u) * (1. - v);
         break;
      case 1:
         if (swap) v *= -1;
         res = 0.5 * deprecatedLegendre::phiFunc(ionedge + 2, v) * (1. + u);
         break;
      case 2:
         if (swap) u *= -1;
         res = 0.5 * deprecatedLegendre::phiFunc(ionedge + 2, u) * (1. + v);
         break;
      case 3:
         if (swap) v *= -1;
         res = 0.5 * deprecatedLegendre::phiFunc(ionedge + 2, v) * (1. - u);
         break;
      case 4:  // Only 3D here...
         //        throw;// Pas code correctement
         if (swap) w *= -1;
         res = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, w) * (1. - u) * (1. - v);
         break;
      case 5:
         if (swap) w *= -1;
         res = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, w) * (1. + u) * (1. - v);
         break;
      case 6:
         if (swap) w *= -1;
         res = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, w) * (1. + u) * (1. + v);
         break;
      case 7:
         if (swap) w *= -1;
         res = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, w) * (1. - u) * (1. + v);
         break;
      case 8:
         if (swap) u *= -1;
         res = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, u) * (1. - v) * (1. + w);
         break;
      case 9:
         if (swap) v *= -1;
         res = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, v) * (1. + u) * (1. + w);
         break;
      case 10:
         if (swap) u *= -1;
         res = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, u) * (1. + v) * (1. + w);
         break;
      case 11:
         if (swap) v *= -1;
         res = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, v) * (1. - u) * (1. + w);
         break;
      default:
         res = 0.;
         break;
   }

   // Correction pour les Hex:
   if (geo_appro->getEntity()->getType() == mEntity::HEX && edge < 4)
   {
      res *= 0.5 * (1. - w);
   }

   return;
}

void xApproxFunctionScalarEdgeHierarchicalLegendreQH_deprecated::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                                         xtensor::xVector<>& res) const
{
   xtensor::xPoint uvw = geo_appro->getUVW();
   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);
   double swapSign = swap ? -1. : 1.;

   // On corrige le cas de l'edge (sinon, on a v et w nuls et donc la fonction de forme aussi)
   if (geo_appro->getEntity()->getType() == mEntity::EDGE)
   {
      v = -1.;
      w = -1.;
   }

   switch (edge)
   {
      case 0:  // 2D et 3D ici...
         if (swap) u *= -1.;
         res[0] = 0.5 * deprecatedLegendre::dPhiFunc(ionedge + 2, u) * (1. - v) * swapSign;
         res[1] = 0.5 * deprecatedLegendre::phiFunc(ionedge + 2, u) * (-1.);
         res[2] = 0;

         if (geo_appro->getEntity()->getType() == mEntity::HEX)
         {
            res[0] *= 0.5 * (1. - w);
            res[1] *= 0.5 * (1. - w);
            res[2] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, u) * (1. - v) * (-1.);
         }
         break;
      case 1:
         if (swap) v *= -1.;
         res[0] = 0.5 * deprecatedLegendre::phiFunc(ionedge + 2, v);
         res[1] = 0.5 * deprecatedLegendre::dPhiFunc(ionedge + 2, v) * (1. + u) * swapSign;
         res[2] = 0;

         if (geo_appro->getEntity()->getType() == mEntity::HEX)
         {
            res[0] *= 0.5 * (1. - w);
            res[1] *= 0.5 * (1. - w);
            res[2] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, v) * (1. + u) * (-1.);
         }
         break;
      case 2:
         if (swap) u *= -1.;
         res[0] = 0.5 * deprecatedLegendre::dPhiFunc(ionedge + 2, u) * (1. + v) * swapSign;
         res[1] = 0.5 * deprecatedLegendre::phiFunc(ionedge + 2, u);
         res[2] = 0;

         if (geo_appro->getEntity()->getType() == mEntity::HEX)
         {
            res[0] *= 0.5 * (1. - w);
            res[1] *= 0.5 * (1. - w);
            res[2] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, u) * (1. + v) * (-1.);
         }
         break;
      case 3:
         if (swap) v *= -1.;
         res[0] = 0.5 * deprecatedLegendre::phiFunc(ionedge + 2, v) * (-1.);
         res[1] = 0.5 * deprecatedLegendre::dPhiFunc(ionedge + 2, v) * (1. - u) * swapSign;
         res[2] = 0.;

         if (geo_appro->getEntity()->getType() == mEntity::HEX)
         {
            res[0] *= 0.5 * (1. - w);
            res[1] *= 0.5 * (1. - w);
            res[2] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, v) * (1. - u) * (-1.);
         }
         break;
      case 4:  // Only 3D here...
         //        throw;
         if (swap) w *= -1.;
         res[0] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, w) * (-1.) * (1. - v);
         res[1] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, w) * (1. - u) * (-1.);
         res[2] = 0.25 * deprecatedLegendre::dPhiFunc(ionedge + 2, w) * (1. - u) * (1. - v);
         break;
      case 5:
         if (swap) w *= -1.;
         res[0] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, w) * (1. - v);
         res[1] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, w) * (1. + u) * (-1.);
         res[2] = 0.25 * deprecatedLegendre::dPhiFunc(ionedge + 2, w) * (1. + u) * (1. - v) * swapSign;
         break;
      case 6:
         if (swap) w *= -1.;
         res[0] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, w) * (1. + v);
         res[1] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, w) * (1. + u);
         res[2] = 0.25 * deprecatedLegendre::dPhiFunc(ionedge + 2, w) * (1. + u) * (1. + v) * swapSign;
         break;
      case 7:
         if (swap) w *= -1.;
         res[0] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, w) * (-1.) * (1. + v);
         res[1] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, w) * (1. - u);
         res[2] = 0.25 * deprecatedLegendre::dPhiFunc(ionedge + 2, w) * (1. - u) * (1. + v) * swapSign;
         break;
      case 8:
         if (swap) u *= -1.;
         res[0] = 0.25 * deprecatedLegendre::dPhiFunc(ionedge + 2, u) * (1. - v) * (1. + w) * swapSign;
         res[1] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, u) * (-1.) * (1. + w);
         res[2] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, u) * (1. - v);
         break;
      case 9:
         if (swap) v *= -1.;
         res[0] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, v) * (1. + w);
         res[1] = 0.25 * deprecatedLegendre::dPhiFunc(ionedge + 2, v) * (1. + u) * (1. + w) * swapSign;
         res[2] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, v) * (1. + u);
         break;
      case 10:
         if (swap) u *= -1.;
         res[0] = 0.25 * deprecatedLegendre::dPhiFunc(ionedge + 2, u) * (1. + v) * (1. + w) * swapSign;
         res[1] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, u) * (1. + w);
         res[2] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, u) * (1. + v);
         break;
      case 11:
         if (swap) v *= -1.;
         res[0] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, v) * (-1.) * (1. + w);
         res[1] = 0.25 * deprecatedLegendre::dPhiFunc(ionedge + 2, v) * (1. - u) * (1. + w) * swapSign;
         res[2] = 0.25 * deprecatedLegendre::phiFunc(ionedge + 2, v) * (1. - u);
         break;
      default:
         res[0] = 0.;
         res[1] = 0.;
         res[2] = 0.;
         break;
   }

   res = (geo_appro)->PushBack(res);
}

void xApproxFunctionScalarQuadHierarchicalLegendre_deprecated::setNodeIndices(int order, int ionface,
                                                                              std::vector<std::vector<int>>& indices,
                                                                              bool tensorSpace)
{
   if (indices.empty())
   {
      std::vector<int> idx(2);

      if (tensorSpace)
      {
         //            cout<<"TENSORSPACE!!!!!\n";
         indices.resize(tensorSpaceDefinition::nbShapeFunctionQuad(order));
         for (int J = 1; J < order; ++J)
         {
            for (int I = 1; I < order; ++I)
            {
               idx[0] = I;
               idx[1] = J;
               indices[(I - 1) + (order - 1) * (J - 1)] = idx;
            }
         }
      }
      else
      {
         //            cout<<"TRUNKSPACE!!!!!\n";
         indices.reserve(trunkSpaceDefinition::nbShapeFunctionQuad(order));
         // Par rapport a Szabo : mes indices (I,J) sont des indices sur le quad, et non pas celui de l'ordre des monomes (i,j)
         // Pour passer de l'un a l'autre, il faut ajouter un.
         for (int J = 1; J < order; ++J)
         {
            int j = J + 1;
            for (int I = 1; I < order; ++I)
            {
               int i = I + 1;
               if (trunkSpaceDefinition::acceptShapeFunc(order, i, j))
               {
                  idx[0] = I;
                  idx[1] = J;
                  indices.push_back(idx);
                  //                        cout<<"pushback "<<i<<" "<<j<<endl;
               }
            }
         }
      }
   }
}

/// applySwap :
/// idx[2] : contient l'indice des variables (AVANT puis APRES swap) : idx[i] = index de la variable i dans uvw
/// signs[2] : contient le signe de la variable : idx[i] = signe de la variable i
/// Exemple : U = -v et V = u :
/// Entree idx=[0,1] car a priori si pas de swap U=u et V=v
/// Sortie idx=[1,0] et signs=[-1,1]

// Precondition : signs initialized to [1,1]
void xApproxFunctionScalarQuadHierarchicalLegendre_deprecated::applySwap(int* idx, double* signs, int swapcase)
{
   int idxIN[2] = {idx[0], idx[1]};

   switch (swapcase)
   {
      case 0:
         // nothing to do
         return;
         break;
      case 1:
         // U = v;
         // V = u;
         idx[0] = idxIN[1];
         idx[1] = idxIN[0];
         //+ signs inchange
         break;
      case 2:
         // U = -u;
         // V = -v;
         signs[0] = -1.;
         signs[1] = -1.;
         //+ idx inchange
         break;
      case 3:
         // U = -u;
         // V = v;
         // idx inchange
         signs[0] = -1.;
         break;
      case 4:
         // U = -v;
         // V = u;
         idx[0] = idxIN[1];
         idx[1] = idxIN[0];
         signs[0] = -1.;
         break;
      case 5:
         // U = -v;
         // V = -u;
         idx[0] = idxIN[1];
         idx[1] = idxIN[0];
         signs[0] = -1.;
         signs[1] = -1.;
         break;
      case 6:
         // U = v;
         // V = -u;
         idx[0] = idxIN[1];
         idx[1] = idxIN[0];
         signs[0] = -1.;
         signs[1] = -1.;
         break;
      case 7:
         // U = u;
         // V = -v;
         // idx inchange
         signs[1] = -1.;
         break;

      default:
         break;
   }

   //    swapDirU = 1;
   //    swapDirV = 1;
}

const int xApproxFunctionScalarQuadHierarchicalLegendre_deprecated::swapvar[8][2] = {{0, 1}, {1, 0}, {1, 0}, {0, 1},
                                                                                     {0, 1}, {1, 0}, {1, 0}, {0, 1}};

xApproxFunctionScalarQuadHierarchicalLegendre_deprecated::xApproxFunctionScalarQuadHierarchicalLegendre_deprecated(
    int _order, int _face, int _ionface, int _swap, bool tensorspace)
    : order(_order), ionface(_ionface), face(_face), swap(_swap)
{
   if (order > XAPPROXFUNCTIONLEGENDREQHMAXORDER)
   {
      cout << "Tabulated functions up to order XAPPROXFUNCTIONLEGENDREQHMAXORDER !\n";
      throw;
   }
   setNodeIndices(order, ionface, indices, tensorspace);
}

void xApproxFunctionScalarQuadHierarchicalLegendre_deprecated::getVal(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                                      double& res) const
{
   //    throw;

   xtensor::xPoint uvw = geo_appro->getUVW();
   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);
   int idx[2] = {0, 1};
   double signs[2] = {1., 1.};

   switch (face)
   {
      case 0:  // 2D et 3D ici...
         // ATTENTION : comme les indices[ionface] commencent a 1 (et pas zero), on a des phifunc(i+1) et pas i+2

         if (swap) applySwap(idx, signs, swap);
         // u, v
         res = deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
               deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]);
         if (geo_appro->getEntity()->getType() == mEntity::HEX) res *= 0.5 * (1. - w);
         //        res = 1.;
         break;
      case 1:  // Only 3D here...
         //        throw;
         idx[1] = 2;
         applySwap(idx, signs, swap);
         //        printf("uu %f u %f ww %f w %f v %f\n",uu,u,ww,w,v);
         // u, w
         res = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
               deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. - v);
         //        res = 2.;
         break;
      case 2:
         idx[0] = 1;
         idx[1] = 2;
         applySwap(idx, signs, swap);
         // v, w
         res = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
               deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. + u);
         //        res = 3.;
         break;
      case 3:
         idx[1] = 2;
         applySwap(idx, signs, swap);
         // u, w
         res = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
               deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. + v);
         //        res = 4.;
         break;
      case 4:
         idx[0] = 1;
         idx[1] = 2;
         applySwap(idx, signs, swap);
         // v, w
         res = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
               deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. - u);
         //        res = 5.;
         break;
      case 5:
         applySwap(idx, signs, swap);
         // u, v
         res = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
               deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. + w);
         //        res = 6.;
         break;
      default:
         res = 0.;
         break;
   }

   return;
}

void xApproxFunctionScalarQuadHierarchicalLegendre_deprecated::getGrad(const xGeomElem* geo_appro, const xGeomElem* geo_integ,
                                                                       xtensor::xVector<>& res) const
{
   xtensor::xPoint uvw = geo_appro->getUVW();
   double u = uvw(0);
   double v = uvw(1);
   double w = uvw(2);

   int idx[2] = {0, 1};
   double signs[2] = {1., 1.};

   switch (face)
   {
      case 0:  // 2D et 3D ici...
         applySwap(idx, signs, swap);
         // u, v
         res[idx[0]] = signs[0] * deprecatedLegendre::dPhiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                       deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]);
         res[idx[1]] = signs[1] * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                       deprecatedLegendre::dPhiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]);
         res[2] = 0.;

         // Correction pour les Hex:
         if (geo_appro->getEntity()->getType() == mEntity::HEX)
         {
            res[idx[0]] *= 0.5 * (1. - w);
            res[idx[1]] *= 0.5 * (1. - w);
            res[2] = -0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                     deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]);
         }
         break;
      case 1:  // Only 3D here...
         //        throw;
         idx[1] = 2;
         applySwap(idx, signs, swap);
         // u, w
         res[idx[0]] = 0.5 * deprecatedLegendre::dPhiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                       deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. - v) * signs[0];
         res[1] = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                  deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (-1.);
         res[idx[1]] = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                       deprecatedLegendre::dPhiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. - v) * signs[1];
         break;
      case 2:
         idx[0] = 1;
         idx[1] = 2;
         applySwap(idx, signs, swap);
         // v, w
         res[0] = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                  deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]);
         res[idx[0]] = 0.5 * deprecatedLegendre::dPhiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                       deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. + u) * signs[0];
         res[idx[1]] = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                       deprecatedLegendre::dPhiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. + u) * signs[1];
         break;
      case 3:
         idx[1] = 2;
         applySwap(idx, signs, swap);
         // u, w
         res[idx[0]] = 0.5 * deprecatedLegendre::dPhiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                       deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. + v) * signs[0];
         res[1] = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                  deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]);
         res[idx[1]] = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                       deprecatedLegendre::dPhiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. + v) * signs[1];
         break;
      case 4:
         idx[0] = 1;
         idx[1] = 2;
         applySwap(idx, signs, swap);
         // v, w
         res[0] = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                  deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (-1.);
         res[idx[0]] = 0.5 * deprecatedLegendre::dPhiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                       deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. - u) * signs[0];
         res[idx[1]] = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                       deprecatedLegendre::dPhiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. - u) * signs[1];
         break;
      case 5:
         applySwap(idx, signs, swap);
         // u, v
         res[idx[0]] = 0.5 * deprecatedLegendre::dPhiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                       deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. + w) * signs[0];
         res[idx[1]] = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                       deprecatedLegendre::dPhiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]) * (1. + w) * signs[1];
         res[2] = 0.5 * deprecatedLegendre::phiFunc(indices[ionface][0] + 1, uvw(idx[0]) * signs[0]) *
                  deprecatedLegendre::phiFunc(indices[ionface][1] + 1, uvw(idx[1]) * signs[1]);
         break;
      default:
         res[0] = 0.;
         res[1] = 0.;
         res[2] = 0.;
         break;
   }

   res = (geo_appro)->PushBack(res);
   return;
}

}  // namespace xfem
