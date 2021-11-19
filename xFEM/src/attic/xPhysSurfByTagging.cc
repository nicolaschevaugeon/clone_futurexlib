/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include <iostream>
#include <sstream>

#include "xPhysSurfByTagging.h"
#include "xPhysSurfParameter.h"

#include "xAttachableChar.h"
#include "xMesh.h"
#include "xLevelSetOperators.h"
#include "xRefCutToAOMD.h"


#ifdef PARALLEL
#include "xParallel.h"
#endif

using AOMD::mEntity;

namespace xfem
{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xcut::xPhysSurfByTagging class implementation ////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////// Constructor ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xcut::xPhysSurfByTagging::xcut::xPhysSurfByTagging(xLevelSet & ls_, xcut::xPhysSurfParameter param )
  : mesh(0),  mesh_bnd(0),region(ls_.getSupport()), classify_in(param.getClassifyerIn()),classify_out(param.getClassifyerOut()), ls(ls_),fittol(param.getFittol()),only_higher_dim_partition(param.getPartitionParameter()), promotors(param.getPromotors())
{

    // setting tags names using level set pointer
    std::ostringstream oss1,oss2;
    oss1 << "tag_entities_xcut::xPhysSurfByTagging_info_"<<&ls_;
    tag_entities = AOMD::AOMD_Util::Instance()->lookupMeshDataId(oss1.str().c_str());
    oss2 << "tag_support_xcut::xPhysSurfByTagging_info_"<<&ls_;
    tag_support = AOMD::AOMD_Util::Instance()->lookupMeshDataId(oss2.str().c_str());

    // init mesh
    mesh = ls.getSupport().getMesh();

    // setting code values
    char *res = (char *) ( &code_support_cover_in );
    res[0] = TE_INCOVER;
    res[1] = SCUTEW_STAT;
    res = (char *) ( &code_support_cover_out );
    res[0] = TE_OUTCOVER;
    res[1] = SCUTEW_STAT;

    // create tagging , cut and classify
    construct_(param.getFitting(), param.getKeepOldPartition(), param.getRecursive());

}
/////////////////////////////////////// End Constructor ////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Destructor /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
xcut::xPhysSurfByTagging::~xcut::xPhysSurfByTagging()
{
    // local variable
    unsigned short int i;
    const int d = ( region.getMesh())->dim()+1;
    xIter itend = region.end(0);
    xIter it = region.end(0);

    // loop on all entities to suppress tagging
    for ( i = 0; i < d; ++i)
    {
        itend = region.end(i);
        for (it = region.begin(i); it != itend; ++it)
        {
            mEntity* e = *it;
            // all entities are tagged with tag_entities so no test on existance of the tag here. We just remove it
            e->deleteData(tag_entities);
            // not perfect, test on pointer would have been more sure, here if a bug on tagging set the value to zero it won't be suppresse
            if (getAttachedChar(e,tag_support))
                e->deleteData(tag_support);
        }
    }

    delete mesh_bnd;


}
/////////////////////////////////////// End Destructor /////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Private methode ////////////////////////////////////////////////////////////////////////////////////////////////////////////
void xcut::xPhysSurfByTagging::construct_(bool fit, bool keep_old_partition_flag, bool recursive)
{

#ifdef TIMING_MONITORING
    xDeltaTime dt;
    int idt = dt.start("level set fitting");
#endif

    if (fit)
    {
        xFitToVertices fit(fittol);
        ls.accept(fit);
    }

#ifdef TIMING_MONITORING
    dt.end(idt);
    idt = dt.start("Cut and tag (seq)");
#endif
    if (mesh_bnd)
        delete mesh_bnd;
    mesh_bnd = new xMesh;

    // generate the iso-zero, cut the mesh and classify sub entity
    // keep_old_partition_flag  indicates that we keep the old partition.
    // first : construction of the object which will do the job
    xcut::xRefCutToAOMD submesh_generator(ls,tag_entities,tag_support,mesh_bnd,classify_in,classify_out,
                                    keep_old_partition_flag,recursive,promotors,only_higher_dim_partition);
    // second : call the methode
    submesh_generator.cutMeshByRef();

#ifdef TIMING_MONITORING
    dt.end(idt);
#endif
    ////////////////////////
    // if in paralle context
#ifdef PARALLEL

#ifdef TIMING_MONITORING
    idt = dt.start("Tag (par)");
#endif

    //  transfer of  tag_entities and tag_support for all type of entity of the frontier of process domaine
    const xPartitionBoundaryInfo &pbi = region.getPartitionBoundaryInfo();
    DataExchange(pbi, *this );

#ifdef TIMING_MONITORING
    dt.end(idt);
#endif

    //  At this stage the exchange have set tempory tag properly according to all process
#endif
    // end in paralle context
    ////////////////////////


    // MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---
    // Here is where we have to put the transfert of tempory tags to public comon tag container
    // It contain all tag of all level-set. See multilevelset  appli for a start. Has now we deal with tag with multiple status,
    // a bit container is no more posible.
    // instead a tab of char is more natural and masking an so one
    // have to be considarate as string comparaison (string base on char tag table). See if possible ...
    // Otherwise bit comparaison with mask of mask ...
    // MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---


#ifdef TIMING_MONITORING
    idt = dt.start("classify");
#endif
    // classify all entity except sub entity using tagging
    // todo :  put this in cutMeshByRef or remove classification in xcut::xPhysSurfByTagging
    int i, d = mesh->dim()+1;
    xIter itend = region.end(0);
    xIter it = region.end(0);
    for ( i = 0; i < d; ++i)
    {
        itend = region.end(i);
        for (it = region.begin(i); it != itend; ++it)
        {
            mEntity* e = *it;
            if (supportCoversIn(e))
                classify_in(e);
            else
                classify_out(e);
        }
    }

#ifdef TIMING_MONITORING
    dt.end(idt);
    dt.print();
#endif

    return;
}

// Promotion methode
void xcut::xPhysSurfByTagging::setPromotorStatus(mEntity* e)
{
   promotor_status=getAttachedChar(e,tag_entities);
}
void xcut::xPhysSurfByTagging::promoteStatus(mEntity* e)
{
   // promoting is mandatory only if entity is not already tagged
   // if it is already tagged it means that it is a sub entity generated during 
   // cutting process of this xcut::xPhysSurf
   if (getAttachedChar(e,tag_entities)==NO_STAT)
          attachChar(e,tag_entities,promotor_status);
}

void xcut::xPhysSurfByTagging::unPromoteStatus(mEntity* e)
{
  e->deleteData(tag_entities);
}

#ifdef PARALLEL

// methode to get and pack tags value of entity e in data container
void xcut::xPhysSurfByTagging::send(mEntity *e, xcut::xPhysSurfByTagging::data_type &data) const
{
    char *res = (char *) ( &data );
    res[0] = getAttachedChar(e,tag_entities);
    res[1] = getAttachedChar(e,tag_support);
    return;

}

// methode to actualise tags value for entity e with information of other process domaine
void xcut::xPhysSurfByTagging::receive(mEntity *e, const std::vector < xcut::xPhysSurfByTagging::data_type > &data) const
{
    char *res;
    char res_cur_entities;
    char res_cur_support;
    bool potentialy_bound = false;
    bool is_rel_out=false;
    bool is_rel_in=false;
    bool potentialy_cut = false;

    res_cur_entities = getAttachedChar(e,tag_entities);
    res_cur_support = getAttachedChar(e,tag_support);


    // here entity tag is normaly consistant betwen process : nothing to do.
    // => only a partial consistant check is done for bug tracking (only when other things have to be done).
    // At some extand the transfert and check sould be removed ...

    // if the entity is allready tagged as "support cut by iso-zero Elt wise" it will stay like that whatever it was tagged in other processors
    // there is nothing to do (no loop)
    if (res_cur_support&SCUTEW_STAT )
        ;
    // if the entity is allready tagged as "support cut by iso-zero " it will stay like that if it was not tagged as "support cut by iso-zero Elt wise"
    // in  other processors
    // => loop needed to check status
    else if (res_cur_support&SCUT_STAT )
    {
        // loop on status of this entity from other process
        std::vector < xcut::xPhysSurfByTagging::data_type >::const_iterator it = data.begin();
        std::vector < xcut::xPhysSurfByTagging::data_type >::const_iterator itend = data.end();
        for (; it != itend; it++)
        {
            res = (char *) ( &( *it ) );

            // check entity status
            if (res[0] != res_cur_entities)
            {
                std::ostringstream oss;
                oss << " Bugg : Inconsitant tagging betwen process, a node already taged as "<<res_cur_entities<<" is tagged "<<res[0]<<" in a other process";
                throw xcut::xRefCutToAOMDException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
            }

            // it was tagged as "cut by iso-zero Elt wise" by an other process
            if (res[1] & SCUTEW_STAT )
            {
                // "cut by iso-zero Elt wise" wins => tag and stop as nothing else have to be donne
                attachChar(e,tag_support,SCUTEW_STAT);
                return;
            }
        }
    }
    // the support have "no" or "other" tagging. Depending on what happend in other process this
    // entity may have a status wich change
    else
    {


        // loop on status of this entity from other process
        std::vector < xcut::xPhysSurfByTagging::data_type >::const_iterator it = data.begin();
        std::vector < xcut::xPhysSurfByTagging::data_type >::const_iterator itend = data.end();
        for (; it != itend; it++)
        {
            res = (char *) ( &( *it ) );

            // check entity status
            if (res[0] != res_cur_entities)
            {
                std::ostringstream oss;
                oss << " Bugg : Inconsitant tagging betwen process, a node already taged as "<<res_cur_entities<<" is tagged "<<res[0]<<" in a other process";
                throw xcut::xRefCutToAOMDException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
            }

            // it was tagged as "cut by iso-zero Elt wise" by an other process
            if (res[1] & SCUTEW_STAT )
            {
                // "cut by iso-zero Elt wise" wins => tag and stop as nothing else have to be donne
                attachChar(e,tag_support,SCUTEW_STAT);
                return;
            }
            // it was tagged as "cut by iso-zero" by an other process
            // Loop must be folowed till the end as "cut by iso-zero  Elt wise"  may appeares for this entity
            else if (res[1] & SCUT_STAT )
            {
                potentialy_cut = true;
            }
            // It was tagged as "in or out touched by iso-zero" by an other process
            //    => it is potentialy in this group
            // Loop must be folowed till the end as "cut by iso-zero  Elt wise"  may appeares for this entity
            else if (res[1] & TS_RELATED )
            {
                potentialy_bound = true;
                if (res[1] &SREL_IN_STAT) is_rel_in=true;
                else is_rel_out=true;
            }


        }

        // after looping on all process contribution it appears that it is realy "cut by iso-zero"
        if (potentialy_cut)
        {
            // "cut by iso-zero " win against "in or out touched by iso-zero" => tag and stop as nothing else have to be donne
            attachChar(e,tag_support,SCUT_STAT);
            return;
        }

        // after looping on all process contribution it appears that it is realy "in or out touched by iso-zero"
        if (potentialy_bound)
        {
            //  check status from this process
            if (res_cur_support&SREL_IN_STAT) is_rel_in=true;
            else if (res_cur_support&SREL_OUT_STAT) is_rel_out=true; 

            //  see if entity is having simultaniously "in touched by iso-zero" and "out touched by iso-zero"  status
            //  => if yes it have to be tagged as "cut by iso-zero" 
            if (is_rel_in && is_rel_out)
               attachChar(e,tag_support,SCUT_STAT);
            else if (is_rel_in)
               attachChar(e,tag_support,SREL_IN_STAT);
            else 
               attachChar(e,tag_support,SREL_OUT_STAT);

        } // end if entity have to be tagged as  "in or out touched by iso-zero"
    } // end else if the entity may have a status changing 

return;

}

#endif

/////////////////////////////////////// End Private methode ////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Public methode /////////////////////////////////////////////////////////////////////////////////////////////////////////////
// update classification and cutting of the domaine (when level set attached to this xcut::xPhysSurf change)
void xcut::xPhysSurfByTagging::update(bool fit, bool keep_old_partition, bool recursive)
{
    construct_(fit, keep_old_partition, recursive);
}

// give dimension of the mesh attached to this xcut::xPhysSurf
int xcut::xPhysSurfByTagging::dim() const {return mesh->dim(); }

// get the mesh or a pointeur to it
xMesh &xcut::xPhysSurfByTagging::getMesh(){return *mesh; }
const xMesh &xcut::xPhysSurfByTagging::getMesh() const {return *mesh; }
xMesh * xcut::xPhysSurfByTagging::getMeshPtr() {return mesh; }
const xMesh * xcut::xPhysSurfByTagging::getMeshPtr() const {return mesh; }

// get the iso-zero mesh
xMesh *  xcut::xPhysSurfByTagging::getMesh_bnd() {return mesh_bnd; }
const xMesh *  xcut::xPhysSurfByTagging::getMesh_bnd() const {return mesh_bnd; }

// Appartenance methode
//
// support
bool xcut::xPhysSurfByTagging::supportCoversIn(mEntity* e) const
{
    unsigned short int code;
    char *res = (char *) ( &code );
    res[0] = getAttachedChar(e,tag_entities);
    res[1] = getAttachedChar(e,tag_support);
    return ( code & code_support_cover_in );
}
bool xcut::xPhysSurfByTagging::supportCoversOut(mEntity* e) const
{
    unsigned short int code;
    char *res = (char *) ( &code );
    res[0] = getAttachedChar(e,tag_entities);
    res[1] = getAttachedChar(e,tag_support);
    return ( code & code_support_cover_out );
}
bool xcut::xPhysSurfByTagging::supportBoundary(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_support)) & TS_BOUNDARY );
}

bool xcut::xPhysSurfByTagging::supportCutStrictly(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_support)) & TS_CUTSTRICTLY );
}

bool xcut::xPhysSurfByTagging::supportCutStrictlyEltWise(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_support)) & SCUTEW_STAT );
}

// atomic
bool xcut::xPhysSurfByTagging::noTouchingIn(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_entities)) & STRICT_IN_STAT );
}
bool xcut::xPhysSurfByTagging::touchingIn(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_entities)) & LOOS_IN_STAT );
}
bool xcut::xPhysSurfByTagging::inIsoZero(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_entities)) & ISO_ZERO_STAT );
}
bool xcut::xPhysSurfByTagging::cutStrictly(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_entities)) & CUT_STAT );
}
bool xcut::xPhysSurfByTagging::touchingOut(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_entities)) & LOOS_OUT_STAT );
}
bool xcut::xPhysSurfByTagging::noTouchingOut(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_entities)) & STRICT_OUT_STAT );
}

// comoposed
bool xcut::xPhysSurfByTagging::strictIn(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_entities)) & TE_INDOMAIN );
}
bool xcut::xPhysSurfByTagging::coversIn(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_entities)) & TE_INDOMAINCUT );
}
bool xcut::xPhysSurfByTagging::strictOut(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_entities)) & TE_OUTDOMAIN );
}
bool xcut::xPhysSurfByTagging::coversOut(mEntity* e) const
{
    return ( ( getAttachedChar(e,tag_entities)) & TE_OUTDOMAINCUT );
}

/////////////////////////////////////// End Public methode /////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xcut::xPhysSurfByTagging class implementation ////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xcut::xPhysSurfByTaggingException class implementation ///////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// general exception used for all xcut::xPhysSurfByTagging throw
xcut::xPhysSurfByTaggingException::xcut::xPhysSurfByTaggingException(std::string info,std::string file,int Line,std::string date,std::string time){
    std::ostringstream oss;
    oss << "In file "<< file << " line " << Line << " compiled "<<date<<" at "<<time<<std::endl;
    oss << "xcut::xPhysSurfByTaggingException : "<< info << std::endl;
    msg = oss.str();
}
/// general exception object : destructor
xcut::xPhysSurfByTaggingException :: ~xcut::xPhysSurfByTaggingException() throw( ) {}

/// mandatory what method
const char * xcut::xPhysSurfByTaggingException::what() const throw( )
{
    return this->msg.c_str();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xcut::xPhysSurfByTaggingException class implementation //////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



} // end of namespace
