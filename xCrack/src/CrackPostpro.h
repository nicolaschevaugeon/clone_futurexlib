/*
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms
    and conditions.
*/

#ifndef _CRACK_postpro_H
#define _CRACK_postpro_H

#include <map>
#include <string>
#include <vector>

// xtensor
#include "xTensor2.h"
#include "xTensor4.h"
#include "xVector.h"

// xfem
#include "xEval.h"
#include "xMesh.h"

// xcrack
#include "lCrack.h"

namespace xcrack
{
xtensor::xTensor2<>& Transpo(xtensor::xTensor2<>& in);

void JGrid3D(const xtensor::xPoint& pointJfront, double rx, double ry, double rz, const xtensor::xTensor2<>& mat,
             std::map<int, xtensor::xPoint>& JCoord, std::map<int, std::vector<int>>& JElements);

void CreateMeshWithPostProElements(xfem::xMesh* mesh, const std::map<int, xtensor::xPoint>& new_coords,
                                   const std::map<int, std::vector<int>>& new_conn3D);

void GetJint3DLevelset(const std::vector<double>& MatInfo, xfem::xEval<xtensor::xTensor2<>>& grad_disp, double alpha,
                       const xtensor::xVector<>& dalpha, const xtensor::xTensor4<>& Phys, double Weight, double DetJac,
                       xfem::xGeomElem* geo_elem, const lCrack* LS_Crack, double& Jhh, xtensor::xVector<>& Ih,
                       xtensor::xTensor2<>& I, double& vol);

void GetJint3DLevelsetThermic(const std::vector<double>& MatInfo, xfem::xEval<xtensor::xTensor2<>>& grad_disp,
                              xfem::xEval<double>& dilatationthermic, double alpha, const xtensor::xVector<>& dalpha,
                              const xtensor::xTensor4<>& Phys, double Weight, double DetJac, xfem::xGeomElem* geo_elem,
                              const lCrack* LS_Crack, double& Jhh, xtensor::xVector<>& Ih, xtensor::xTensor2<>& I, double& vol);

/// return for a given point values off the assymptotic field for on of the opening mode
/*! on Entries :
    Xloc : Position off the point in a frame link to the crack. The 0 is on the crack tip.
           the Z axis is the local crack lip direction.
           X and Y define a plan orthogonal to the crack plane.
           X is in the crack plan, pointing outside the crack.
    MatInfo : vector of double containing Young Modulus and Poisson ration.
              MatInfo[0] = E
              MatInfo[1] = nu
    mode : openning mode to do the computation should be 1, 2 or 3.

    Output :
    AuxStress : Cauchy Stress Tensor.
    AuxEps    : Symetric despacement Gradient
    AuxDips   : Displacement Vector
    AuxGrad\i\Phys : Gradiant in direction i (refering to X and Y of the above) off physical quantity Phys
    AuxGrad\i\j\Disp : Second order derivative i and j of the Displacement.

!*/
void GetMech3DAuxFields(const xtensor::xPoint& Xloc, const std::vector<double>& MatInfo, xtensor::xTensor2<>& AuxStress,
                        xtensor::xTensor2<>& AuxEps, xtensor::xVector<>& AuxDisp, xtensor::xVector<>& AuxGrad0Disp,
                        xtensor::xVector<>& AuxGrad1Disp, xtensor::xTensor2<>& AuxGrad0Stress,
                        xtensor::xTensor2<>& AuxGrad1Stress, xtensor::xVector<>& AuxGrad00Disp, xtensor::xVector<>& AuxGrad11Disp,
                        xtensor::xVector<>& AuxGrad01Disp, int mode);

void GetGradqLoc(int postid, const xtensor::xPoint& xyz, xtensor::xTensor2<>& GradqLoc, xtensor::xPoint& xyzloc,
                 const std::map<int, xtensor::xVector<>>& JBoxDimension, const std::map<int, xtensor::xPoint>& JFrontPoint,
                 const std::map<int, xtensor::xTensor2<>>& JMap);

double GetqLoc(int postid, const xtensor::xPoint& xyz, const std::map<int, xtensor::xVector<>>& JBoxDimension,
               const std::map<int, xtensor::xPoint>& JFrontPoint, const std::map<int, xtensor::xTensor2<>>& JMap);

double AngleWithMaxHoopStress(const double& K1, const double& K2);

void CreateMeshWithPostProElements(xfem::xMesh* mesh, const std::map<int, xtensor::xPoint>& new_coords,
                                   const std::map<int, std::vector<int>>& new_conn3D);

}  // namespace xcrack

#endif
