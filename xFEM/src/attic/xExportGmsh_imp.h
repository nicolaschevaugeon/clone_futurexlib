/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef _export_gmsh_imp_H
#define _export_gmsh_imp_H


template <typename O>
O & operator << (O & ofs, const xGmshDrawable &p)
{
  switch(p.type)
    {
    case xGmshDrawable::NONE:return ofs;
    case xGmshDrawable::SP  :ofs << "SP (";break;
    case xGmshDrawable::SL  :ofs << "SL (";break;
    case xGmshDrawable::ST  :ofs << "ST (";break;
    case xGmshDrawable::SS  :ofs << "SS (";break;
    case xGmshDrawable::SQ  :ofs << "SQ (";break;
    case xGmshDrawable::SH  :ofs << "SH (";break;
    case xGmshDrawable::VP  :ofs << "VP (";break;
    case xGmshDrawable::VL  :ofs << "VL (";break;
    case xGmshDrawable::VT  :ofs << "VT (";break;
    case xGmshDrawable::VS  :ofs << "VS (";break;
    case xGmshDrawable::VQ  :ofs << "VQ (";break;
    case xGmshDrawable::VH  :ofs << "VH (";break;
    case xGmshDrawable::TP  :ofs << "TL (";break;
    case xGmshDrawable::TL  :ofs << "TL (";break;
    case xGmshDrawable::TT  :ofs << "TT (";break;
    case xGmshDrawable::TS  :ofs << "TS (";break;
    case xGmshDrawable::TQ  :ofs << "TQ (";break;
    case xGmshDrawable::TH  :ofs << "TH (";break;
    default: assert(0); break;
    }


  for(int i=0;i<p.nbNod;i++)
    {
      if(i)ofs << "," << p.X[i] << "," << p.Y[i] << "," << p.Z[i];
      else ofs << p.X[i] << "," << p.Y[i] << "," << p.Z[i];
    }
  ofs << "){";
  for(int i=0;i<p.nbNod*p.nbSamples*p.sizeInfo;i++)
    {
      if(i) ofs << "," << p.V[i];
      else ofs << p.V[i];
    }
  ofs << "};\n";
  
  return ofs;
}

#endif
