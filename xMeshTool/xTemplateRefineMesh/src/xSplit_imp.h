/*
   This file is a part of eXlibris C++ Library
   under the GNU Lesser General Public License.
   See the NOTICE & LICENSE files for terms
   and conditions.
 */

#ifndef XSPLIT_IMP_H
#define XSPLIT_IMP_H

namespace xmeshtool
{
//===== xMeshSplitUserCriteria ===================================================================================================
template < typename T, template < typename > class DM >
xMeshSplitUserCriteria < T,DM >::xMeshSplitUserCriteria(criteria_t &criteria_, int max_it_split_ ) :
    warper(nullptr)
    ,criteria(criteria_)
    ,target_dim(0)
    ,max_it_split(max_it_split_)
    ,part_man(nullptr)
    ,world(nullptr)
    ,new_hanging_key_cont(nullptr)
    ,kmso(nullptr)
    ,kmsar(nullptr)
{
    treat_new_hanging = std::bind( &xMeshSplitUserCriteria < T,DM >::hangingSynchro,std::ref(*this),std::placeholders::_1);
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
void xMeshSplitUserCriteria < T,DM >::registerTransfer(transfert_info_t* trans_information)
{
    trans_information_container.push_front(trans_information);
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
void xMeshSplitUserCriteria < T,DM >::clearInternal()
{
    hanging_down.clear();
    hanging_up.clear();
    trans_information_container.clear();
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
bool xMeshSplitUserCriteria < T,DM >::split(warper_t &warp,int target_dim_, partman_t & part_man_)
{
    //set information used in private methode
    target_dim = target_dim_;
    warper = &warp;
    part_man = &part_man_;
    MPI_Comm world_tmp = part_man->getComm();
    world = &world_tmp;
    int proc_id;
    MPI_Comm_rank(world_tmp,&proc_id);
    switch (target_dim)
    {
        case 1 :
         {
             select_hanging_edge = std::bind( &xMeshSplitUserCriteria < T,DM >::dontSelect,std::ref(*this),std::placeholders::_1);
             select_hanging_face = std::bind( &xMeshSplitUserCriteria < T,DM >::dontSelect,std::ref(*this),std::placeholders::_1);
             delete_hanging_edge_info_from_tet = std::bind( &xMeshSplitUserCriteria < T,DM >::dontDeleteEntityWithCondition,std::ref(*this),std::placeholders::_1,std::placeholders::_2);
             delete_hanging_edge_info_from_face = std::bind( &xMeshSplitUserCriteria < T,DM >::dontDeleteEntityWithCondition,std::ref(*this),std::placeholders::_1,std::placeholders::_2);
             delete_hanging_face_info = std::bind( &xMeshSplitUserCriteria < T,DM >::dontDeleteEntityWithCondition,std::ref(*this),std::placeholders::_1,std::placeholders::_2);
             split_elem = std::bind( &xMeshSplitUserCriteria < T,DM >::splitEdge,std::ref(*this),std::placeholders::_1);
             //split_related = std::bind( &xMeshSplitUserCriteria < T,DM >::xxxxxxxxx,std::ref(*this),std::placeholders::_1);
             break;
         }
        case 2 :
         {
             select_hanging_edge = std::bind( &xMeshSplitUserCriteria < T,DM >::doSelect,std::ref(*this),std::placeholders::_1);
             select_hanging_face = std::bind( &xMeshSplitUserCriteria < T,DM >::dontSelect,std::ref(*this),std::placeholders::_1);
             delete_hanging_edge_info_from_tet = std::bind( &xMeshSplitUserCriteria < T,DM >::dontDeleteEntityWithCondition,std::ref(*this),std::placeholders::_1,std::placeholders::_2);
             delete_hanging_edge_info_from_face = std::bind( &xMeshSplitUserCriteria < T,DM >::deleteEdgeWithCondition,std::ref(*this),std::placeholders::_1,std::placeholders::_2);
             delete_hanging_face_info = std::bind( &xMeshSplitUserCriteria < T,DM >::dontDeleteEntityWithCondition,std::ref(*this),std::placeholders::_1,std::placeholders::_2);
             split_elem = std::bind( &xMeshSplitUserCriteria < T,DM >::splitFace,std::ref(*this),std::placeholders::_1);
             //split_related = std::bind( &xMeshSplitUserCriteria < T,DM >::xxxxxxxxx,std::ref(*this),std::placeholders::_1);
             break;
         }
        case 3 :
         {
             select_hanging_edge = std::bind( &xMeshSplitUserCriteria < T,DM >::doSelect,std::ref(*this),std::placeholders::_1);
             select_hanging_face = std::bind( &xMeshSplitUserCriteria < T,DM >::doSelect,std::ref(*this),std::placeholders::_1);
             delete_hanging_edge_info_from_tet = std::bind( &xMeshSplitUserCriteria < T,DM >::deleteEdgeWithCondition,std::ref(*this),std::placeholders::_1,std::placeholders::_2);
             delete_hanging_edge_info_from_face = std::bind( &xMeshSplitUserCriteria < T,DM >::dontDeleteEntityWithCondition,std::ref(*this),std::placeholders::_1,std::placeholders::_2);
             delete_hanging_face_info = std::bind( &xMeshSplitUserCriteria < T,DM >::deleteFaceWithCondition,std::ref(*this),std::placeholders::_1,std::placeholders::_2);
             split_elem = std::bind( &xMeshSplitUserCriteria < T,DM >::splitTet,std::ref(*this),std::placeholders::_1);
             split_related = std::bind( &xMeshSplitUserCriteria < T,DM >::splitTetRelatedToHanging,std::ref(*this),std::placeholders::_1);
             break;
         }
    }

    // local
    int iter = 1;
    int nbcut_tot,bnd_f;
    const bool verbose = true;
    std::set < entity_t * > to_split_onelevel;
    std::set < entity_t * > bnd_forced_to_split;
    std::set < entity_t * > bnd_forced_to_split_treated;
    std::set < entity_t * > bnd_which_add_to_to_split;

    // Key manager
    kmsar = new splitKeyManagerSendAndReceive < partman_t, entity_t >(*part_man);
    kmso = new splitKeyManagerSendOnly < partman_t, entity_t >(*part_man);

    // Key container
    xtool::xKeyContainerSendOrRecv < information_key_t > one_level_key_cont(*world);
    new_hanging_key_cont = new xtool::xKeyContainerSendOrRecv < information_key_t >(*world);
    xtool::xKeyContainerSendAndRecv < information_key_t > hanging_update_PM_FE_key_cont(*world);

    // Information manager
    splitInfoManagerOneLevelInformation < partman_t, entity_t > one_level_info( *part_man,bnd_which_add_to_to_split);
    splitInfoManagerNewHangingTreatement < partman_t, entity_t > new_hanging_info(part_man_,treat_new_hanging);
    splitInfoManagerUpdatePMFE < partman_t, entity_t, DM < entity_t * >, DM < t1_t >, DM < t2_t > > hanging_update_PM_FE_info( part_man_,hanging_down,hanging_up,hanging_data_t1,hanging_data_t2);

    // global loop on split pass
    do
    {
        if (verbose)
            std::cout << "Entering spliting iteration "<<iter<<std::endl;
        nbcut_tot = 0;
        range_t range_elem(warper->begin(target_dim),warper->end(target_dim));


        // update criteria against new mesh
        criteria.update();

        // reset new hanging keys
        new_hanging_key_cont->clearKeys();


        // create list of element to split in this pass
        // NOTE : except for messages to_split is allways used in way compatible with forward list. Switching
        // to this container would be kool
        std::list < entity_t * > to_split;
        for (auto e : range_elem)
            if (criteria(*e)) to_split.push_front(e);

        // loop till list of element to split is fully completed
        // by one-level rule on boundary
        do
        {
            if (verbose)
            {
                std::cout << "Entering one level loop\nTo be split :  "<<to_split.size()<<std::endl;
            }
            //reset to null
            bnd_forced_to_split.clear();

            // Loop till list of element to split is fully completed by one-level rule.
            // In this loop bnd element (connecting 2 or more elements in between procs) are collected
            // when related local element impose one-level rule to remotes.
            do
            {
                // check one level criteria
                to_split_onelevel.clear();
                for (auto e : to_split) checkOneLevelCriteria(e,to_split_onelevel,bnd_forced_to_split);

                // clear already present entity
                for (auto e : to_split) to_split_onelevel.erase(e);

                // Add entity added by one level criteria into to_split, first.
                // Order first here is a key algorithmic point. It insure that element are treated in
                // reverse order of dependence : i.e. first split element are not concerned by criteria
                // but there treatment permit to split those who are concerned without breaking the one-level rule.
                to_split.insert(to_split.begin(),to_split_onelevel.begin(),to_split_onelevel.end());

            } while (to_split_onelevel.size());

            if (verbose)
            {
                std::cout << "After adding local 1-level impacted element :  "<<to_split.size()<<std::endl;
            }

            // clear already treated bnd
            // This is to avoid unnecessary communication
            // Why is it kept for all level honestly I don't remenber. TODO check this algorithmic aspect.
            for (auto e : bnd_forced_to_split_treated) bnd_forced_to_split.erase(e);

            // communicate on all  boundary to send potential 1 level problem to each other
            bnd_which_add_to_to_split.clear();
            one_level_key_cont.clearKeys();
            one_level_key_cont.accumulateKeys( bnd_forced_to_split.begin(), bnd_forced_to_split.end(),*kmso );
            exchangeInformation(one_level_key_cont,one_level_info);

            // add bnd exchanged to treated ones
            bnd_forced_to_split_treated.insert(bnd_forced_to_split.begin(),bnd_forced_to_split.end());

            // add new new 1 level split to to_split
            addOneLevelCriteria(bnd_which_add_to_to_split,to_split);
            if (verbose)
            {
                std::cout << "After adding  remote 1-level impacted element :  "<<to_split.size()<<std::endl;
            }


            // Check if there is still something to treat. If no bnd element are received from other proc all
            // local one-level rule been applied and remote does not impose anything, it's finish
            // If at least one bnd element is received new related element have to be split and this
            // may impose new extra one-level rule splitting : loop
            bnd_f = bnd_forced_to_split.size();
            // get maximum to loop the same way on all proc (exchange requirement : every body are participating but only
            // some are communicating in fact)
            MPI_Allreduce(MPI_IN_PLACE,&bnd_f,1,MPI_INT,MPI_MAX,*world);

        } while (bnd_f);

        // split entity
        for (auto pe : to_split)
            if (split_elem(*pe)) ++nbcut_tot;

        // get maximum to loop the same way on all proc (exchange mandatory)
        MPI_Allreduce(MPI_IN_PLACE,&nbcut_tot,1,MPI_INT,MPI_MAX,*world);

        // DEBUG DEBUG
        //if (verbose)
        //{
        //    //warper->printMesh(true);
        //    warper->writeToFile("imesh_level_"+std::to_string(iter)+"_proc_"+std::to_string(proc_id)+".msh");
        //}
        // DEBUG DEBUG

        // Exchange new created hanging  information and create on both sides new hanging face/edge counterpart.
        // This feed edge_face_to_be_treated so it has to be reset.
        edge_face_to_be_treated.clear();
        exchangeInformation(*new_hanging_key_cont,new_hanging_info);

        hanging_update_PM_FE_key_cont.clearKeys();
        hanging_update_PM_FE_key_cont.accumulateKeysAllGather( edge_face_to_be_treated.begin(), edge_face_to_be_treated.end(), *kmsar);
        exchangeInformation(hanging_update_PM_FE_key_cont,hanging_update_PM_FE_info);

        // treatement of dangling edges/faces
        removeDangling();

#ifndef NDEBUG
        // save intermediate mesh
        if (verbose)
        {
            warper->writeToFile("mesh_level_"+std::to_string(iter)+"_proc_"+std::to_string(proc_id)+".msh");
        }
#endif
    } while (++iter < max_it_split && nbcut_tot);

    // reset new hanging keys
    new_hanging_key_cont->clearKeys();
    edge_face_to_be_treated.clear();

    // treat now terminal element : element that have a edge/face supporting hanging entity
    // do a copy has warper db will be modyfied
    std::vector < entity_t * > elem_to_treat(warper->begin(target_dim),warper->end(target_dim));
    for (auto e : elem_to_treat)
    {
        // here returned value may be used to store e  vertex to be able to clean
        // up latter terminal element and start again splitting
        split_related(*e);
    }
    // note : split_related in 3D split faces (when all its surounded edges are split). This
    // imply that new hanging entity may have to be treated. Only faces are concerned

    // Exchange new created hanging  information and create on both sides new hanging face/edge counterpart.
    // update partition manager for new hanging edges/faces
    hanging_update_PM_FE_key_cont.clearKeys();
    hanging_update_PM_FE_key_cont.accumulateKeysAllGather( edge_face_to_be_treated.begin(), edge_face_to_be_treated.end(), *kmsar);
    exchangeInformation(hanging_update_PM_FE_key_cont,hanging_update_PM_FE_info);

    // treatement of dangling edges/faces
    removeDangling();

    // cleaning
    if (new_hanging_key_cont ) delete new_hanging_key_cont;
    if (kmso ) delete kmso;
    if (kmsar ) delete kmsar;
    new_hanging_key_cont = nullptr;
    kmsar = nullptr;
    kmso = nullptr;

    return true;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
void xMeshSplitUserCriteria < T,DM >::checkOneLevelCriteria( entity_t *e, std::set < entity_t * > &to_split_onelevel, std::set < entity_t * > &bnd_forced_to_split)
{
    // loop on sub entity to see which are hanging
    for (int i = 1; i < target_dim; ++i)
    {
        assert(e->size(i));
        for (int j = 0, n = e->size(i); j < n; ++j)
        {
            entity_t *se = (entity_t *) e->get(i,j);
            assert(se);

            // try to get upper entity
            entity_t **phe = hanging_up.getData(*se);
            // if there is a uper entity this edge/face is already hanging with a
            // edge/face. Cutting element e will introduce 2 level of hanging.
            // Neighbors have to be split to remove hanging situation so that element e may be split
            if ( phe )
            {
                entity_t *he = *phe;
                // If hanging entity is on proc boundary its remote neighbor have to be split
                // store this edge/face to communicate with others latter
                if (part_man->hasRemoteCopy(*he))
                    bnd_forced_to_split.insert(he);
                // If  hanging entity is local just had neighbor element to to_split_onelevel
                else
                    to_split_onelevel.insert(he->begin(target_dim),he->end(target_dim));
            }
        }
    }
    return;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
void xMeshSplitUserCriteria < T,DM >::addOneLevelCriteria(std::set < entity_t * > &bnd_which_add_to_to_split, std::list < entity_t * > &to_split)
{
    // add upper entity related to add list
    std::set < entity_t * > tmp;
    for (auto se : bnd_which_add_to_to_split)
    {
        // if things are correctly done this entity must already have a associated downward hanging
        // entity
        assert( hanging_down.getData(*se) );

        // loop to store neighbors element of se into tmp
        // se->size(target_dim) may be null if in one proc it is dangling and in an other it is not
        // in the dangling case it is null
        for (int l = 0; l < se->size(target_dim); ++l)
        {
            entity_t *up = se->get(target_dim,l);
            tmp.insert(up);
        }
    }

    // remove already chosen element to split
    for (auto e : to_split)
    {
        tmp.erase(e);
    }

    // add in front new ones
    to_split.insert(to_split.begin(),tmp.begin(),tmp.end());

    return;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
auto xMeshSplitUserCriteria < T,DM >::splitTet(entity_t &t)->entity_t *
{
    // top entity. It may not be considered hanging on something bigger

    //local
    int i;
    // local copies of edges and vertex
    entity_t *e[6];
    entity_t *v[4];

    // checking/filling
    for (i = 0; i < 6; ++i)
    {
        e[i] = t.get(1,i);
        // if this edge is hanging (i.e. is on a biger edge related to something) then
        // something went wrong in onelevel selection. We pospone this tet split
        if ( hanging_up.getData(*e[i]) )
            return nullptr;
        if (  i < 4)
        {
            // retrive tet nodes
            v[i] = t.get(0,i);

            // if this face is hanging (i.e. is on a biger face related to something) then
            // something went wrong in onelevel selection. We pospone this tet split
            if ( hanging_up.getData(*( t.get(2,i))) )
                return nullptr;
        }
    }

    // Treat faces
    entity_t *subface[4];
    std::vector < entity_t * > faces(t.begin(2),t.end(2));    // has we delete adjancy element in the loop we need a copy of adjancy to be correct
    for ( i = 0; i < 4; ++i)
    {
        entity_t &f = *faces[i];
        // split face and get midle sub face
        entity_t *fc = splitFace(f);
        assert (fc);
        subface[i] = fc;

        // remove hanging information if possible
        delete_hanging_face_info(f,*fc);
    }
    // Retrieve mid nodes and remove dangling edge information if not on bnd
    entity_t *ve[6];
    for (i = 0; i < 6; ++i)
    {
        entity_t **pve = hanging_down.getData(*e[i]);
        assert(pve);
        ve[i] = *pve;

        // remove hanging information if possible
        // if not dangling edge treat connectivity TODO ??? why removeEntityWithAdj do the job normally ???
        if ( !delete_hanging_edge_info_from_tet(*e[i],*ve[i]) )
        {
            e[i]->del(&t);
        }
    }

    // collect information related to this tet (not done befor has splitFace use this same mechanism)
    collect(&t);

    // Table for loop describing base face of corner tet
    const int i1[4] = {0,1,2,3};
    const int i2[4] = {2,0,1,5};
    const int i3[4] = {3,4,5,4};

    // Table for loop describing node opposite to base face of corner tet for each internal tet : 3 set coresponding
    // respectively to 3 possible choice of diamond internal egde (3-1,5-0,2-4)
    const int i8[4] = {1,3,3,1};
    const int i9[4] = {5,5,0,0};
    const int i10[4] = {4,2,4,2};

    // Loop to compute max aspect ratio for each internal edge solution
    double ar,ar31 = 0.,ar50 = 0.,ar24 = 0.;
    for (i = 0; i < 4; ++i)
    {
        ar = warper->tetAspectRatio(ve[i1[i]],ve[i2[i]],ve[i3[i]],ve[i8[i]]);
        ar31 = ( ar31 > ar ) ? ar31 : ar;
        ar = warper->tetAspectRatio(ve[i1[i]],ve[i2[i]],ve[i3[i]],ve[i9[i]]);
        ar50 = ( ar50 > ar ) ? ar50 : ar;
        ar = warper->tetAspectRatio(ve[i1[i]],ve[i2[i]],ve[i3[i]],ve[i10[i]]);
        ar24 = ( ar24 > ar ) ? ar24 : ar;
    }


    // find smallest aspect ratio among all choice
    char choice = ( ar31 > ar50 ) ? (( ar50 < ar24 ) ? 2 : 3 ) : (( ar31 < ar24 ) ? 1 : 3 );

    // Table for loop to chose edges and faces
    int ii0,ii1;
    int i4[4], i5[4], i6[4], i7[4];

    //choice=1;

    // set  data according to diamon edge chosen
    switch (choice)
    {
        case 1 :
         {
             ii0 = 3;
             ii1 = 1;
             i4[0] = 0; i4[1] = 4; i4[2] = 5; i4[3] = 2;
             i5[0] = 0; i5[1] = 1; i5[2] = 3; i5[3] = 2;
             i6[0] = 0; i6[1] = 0; i6[2] = 2; i6[3] = 2;
             i7[0] = 3; i7[1] = 1; i7[2] = 3; i7[3] = 1;
             break;
         }
        case 2 :
         {
             ii0 = 0;
             ii1 = 5;
             i4[0] = 3; i4[1] = 2; i4[2] = 1; i4[3] = 4;
             i5[0] = 3; i5[1] = 2; i5[2] = 0; i5[3] = 1;
             i6[0] = 1; i6[1] = 2; i6[2] = 1; i6[3] = 3;
             i7[0] = 0; i7[1] = 3; i7[2] = 2; i7[3] = 0;
             break;
         }
        case 3 :
         {
             ii0 = 2;
             ii1 = 4;
             i4[0] = 3; i4[1] = 5; i4[2] = 1; i4[3] = 0;
             i5[0] = 1; i5[1] = 0; i5[2] = 2; i5[3] = 3;
             i6[0] = 0; i6[1] = 2; i6[2] = 2; i6[3] = 1;
             i7[0] = 3; i7[1] = 3; i7[2] = 1; i7[3] = 0;
             break;
         }
    }

    // create internal edge chosen
    entity_t *ei = warper->createEdge(ve[ii0],ve[ii1]);
    transfer(ei);       // ei inherit original tet property has it is inside it.

    // four inside face linked with internal edge
    entity_t *ifs[4];
    for (i = 0; i < 4; ++i)
    {
        //inside face
        entity_t *e3 = warper->getEdgeFromVertex(ve[ii1],ve[i4[i]]);
        assert(e3);
        entity_t *e2 = warper->getEdgeFromVertex(ve[ii0],ve[i4[i]]);
        assert(e2);
        ifs[i] = warper->createFaceWithEdge(ei,e2,e3);
        transfer(ifs[i]);   // ifs inherit original tet property has it is inside it.
    }
    // tet creation
    entity_t *nt;
    for (i = 0; i < 4; ++i)
    {
        //inside face for corner tet
        entity_t *e1 = warper->getEdgeFromVertex(ve[i1[i]],ve[i2[i]]);
        assert(e1);
        entity_t *e2 = warper->getEdgeFromVertex(ve[i2[i]],ve[i3[i]]);
        assert(e2);
        entity_t *e3 = warper->getEdgeFromVertex(ve[i3[i]],ve[i1[i]]);
        assert(e3);
        entity_t *f2 = warper->createFaceWithEdge(e1,e2,e3);
        transfer(f2);   // f2 inherit original tet property has it is inside it.
        // Lateral faces
        entity_t *f0 = warper->getFaceFromVertex(ve[i2[i]],v[i],ve[i1[i]]);
        entity_t *f1 = warper->getFaceFromVertex(ve[i3[i]],v[i],ve[i1[i]]);
        entity_t *f3 = warper->getFaceFromVertex(ve[i3[i]],v[i],ve[i2[i]]);
        // corner tet
        nt = warper->createTetWithFace(f0,f1,f2,f3);
        transfer(nt);
        // tet connected with internal edge and corner tet
        nt = warper->createTetWithFace(subface[i5[i]],f2,ifs[i6[i]],ifs[i7[i]]);
        transfer(nt);
    }

    warper->removeEntityWithAdj(t);
    return nt;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
auto xMeshSplitUserCriteria < T,DM >::splitFace(entity_t &face)->entity_t *
{

    // if face already split return middle hanging face
    entity_t **pfc = hanging_down.getData(face);
    if (pfc)
        return *pfc;

    // face here are not supose to be anything else then triangle (i.e quadrangle not supported)
    assert(face.size(1) == 3);

    // split edges and store central vertex
    entity_t *ve[4];
    for (int i = 0; i < 3; ++i)
    {
        entity_t &e = *face.get(1,i);
        ve[i] = splitEdge(e);
        delete_hanging_edge_info_from_face(e,*ve[i]);
    }

    // collect information related to this face (not done befor has splitEdge use this same mechanism)
    collect(&face);

    // new internal edges
    entity_t *newe[3];
    for (int i = 0; i < 3; ++i)
    {
        const int i1 = ( i+1 )%3;
        newe[i] = warper->createEdge(ve[i],ve[i1]);
        hanging_up.setData(*newe[i]) = &face;
        transfer(newe[i]);
    }
    // new central face
    entity_t *fc = warper->createFaceWithEdge(newe[0],newe[1],newe[2]);
    hanging_down.setData(face) = fc;
    hanging_up.setData(*fc) = &face;
    transfer(fc);

    // new faces around central one
    for (int i = 0; i < 3; ++i)
    {
        const int i1 = ( i+1 )%3;
        entity_t *vc = face.get(0,i1);
        entity_t *e1 = warper->getEdgeFromVertex(vc,ve[i]);
        assert(e1);
        entity_t *e2 = warper->getEdgeFromVertex(vc,ve[i1]);
        assert(e2);
        entity_t *f = warper->createFaceWithEdge(e1,e2,newe[i]);
        hanging_up.setData(*f) = &face;
        transfer(f);
    }

    // select for hanging treatement
    select_hanging_face(&face);

    return fc;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
auto xMeshSplitUserCriteria < T,DM >::splitEdge(entity_t &edge)->entity_t *
{

    // if edge already split return middle hanging vertex
    entity_t **pv = hanging_down.getData(edge);
    if (pv)
        return *pv;


    // edge vetex
    entity_t *v1 = edge.get(0,0);
    entity_t *v2 = edge.get(0,1);
    // new vertex in the middle of the edge
    entity_t *newv = warper->createCOGVertex(edge);

    // set has hanging down for edge
    hanging_down.setData(edge) = newv;

    // new edges
    // NOTE : order is used in other places, don't change without checking
    entity_t *nedj1 = warper->createEdge(v1,newv);
    entity_t *nedj2 = warper->createEdge(v2,newv);

    // set has hanging up for new edges
    hanging_up.setData(*nedj1) = &edge;
    hanging_up.setData(*nedj2) = &edge;

    // transfer information
    collect(&edge);
    transfer(newv);
    transfer(nedj1);
    transfer(nedj2);

    // select for hanging treatement
    select_hanging_edge(&edge);

    return newv;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
auto xMeshSplitUserCriteria < T,DM >::splitTetRelatedToHanging(entity_t  &t)->entity_t *
{
    int i,j,k,idc = 0;

    // test if hanging edge exist
    // note : hanging face imply hanging edge exist => no test on face
    assert(t.size(1) == 6);
    assert(t.size(0) == 4);
    // uncoment to create a COG of all nodes used for tet (i.e. nodes of original tet + nodes of splited edges)
    //warper->resetPointCoord();
    for (i = 0, k = 0; i < 6; ++i)
    {
        entity_t **pv = hanging_down.getData(*t.get(1,i));
        if ( pv )
        {
            k++;
            idc = i;
            // uncoment to create a COG of all nodes used for tet (i.e. nodes of original tet + nodes of splited edges)
            //warper->accumulatePointCoord(*pv);
        }
    }
    // no central node by default
    entity_t *newv = nullptr;

    // hanging edge exist
    if (k)
    {
        bool trans[10] = {0,0,0,0,0,0,0,0,0,0};
        entity_t *tet[16];
        entity_t *v[4];

        //std::cout<<" hanging: tet ini "<<warper->tetAspectRatio(&t);
        switch (k)
        {
            // only one edge is split
            //case 1 :
            case 99999 : // for now only the central node strategy have to be used.
                         // This part is implemented and run apparently well but no serious
                         // checking have been done. Uncomment to make it work but it is at your risk
                         // It offers better element shape and less elements => to be preferred when checked
                         // But this dos not solve the aspect ratio problem for k=2,3, ... Those have to be
                         // also implemented in a nice solution way (no central node ??)
                         //
             {
                 const int ef[6][13] = {                                                  //idc
                     //  0   1   2    3   4   5    6   7    8    9    10   11  12
                     // f1  f2  vs1 vs2  v0  v1  edg2 edg1 edg3 edg4 edg5 fa2 fa3
                     {0,  1,  2,   3,  0,  1,   2,   1,   3,   4,  5,   2,  3  },        //0
                     {0,  2,  0,   3,  1,  2,   0,   2,   4,   5,  3,   3,  1  },        //1
                     {0,  3,  1,   3,  2,  0,   1,   0,   5,   3,  4,   1,  2  },        //2
                     {1,  3,  1,   2,  0,  3,   0,   4,   2,   5,  1,   2,  0  },        //3
                     {1,  2,  0,   2,  3,  1,   3,   0,   5,   1,  2,   0,  3  },        //4
                     {2,  3,  1,   0,  3,  2,   4,   1,   3,   2,  0,   0,  1  }         //5
                 };
                 //std::cout<<" idc is : "<<idc<<std::endl;
                 warper->printEntity(&t,true);
                 entity_t **pv = hanging_down.getData(*t.get(1,idc));
                 entity_t *nvc = *pv;
                 entity_t *vs1 = t.get(0,ef[idc][2]);
                 entity_t *vs2 = t.get(0,ef[idc][3]);
                 entity_t *v0 = t.get(0,ef[idc][4]);
                 entity_t *v1 = t.get(0,ef[idc][5]);

                 // edges cut
                 assert(warper->getEdgeFromVertex(nvc,v0));
                 assert(warper->getEdgeFromVertex(nvc,v1));
                 entity_t *ec0 = warper->getEdgeFromVertex(nvc,v0);
                 entity_t *ec1 = warper->getEdgeFromVertex(nvc,v1);

                 warper->printEntity(vs1,false);
                 warper->printEntity(vs2,false);
                 warper->printEntity(v0,false);
                 warper->printEntity(v1,false);
                 warper->printEntity(ec0,false);
                 warper->printEntity(ec1,false);

                 //====== First face treatement
                 entity_t *f1 = t.get(2,ef[idc][0]);
                 assert(f1->size(0) == 3);

                 // collect information related to first face
                 collect(f1);

                 // create new edge
                 entity_t *e1 = warper->getEdgeFromVertex(nvc,vs1);
                 if (!e1)
                 {
                     e1 = warper->createEdge(nvc,vs1);
                     transfer(e1);
                 }
                 warper->printEntity(e1,false);
                 // create new face
                 entity_t *f10 = warper->getFaceFromVertex(nvc,vs1,v0);
                 if (!f10)
                 {

                     f10 = warper->createFaceWithEdge(e1,t.get(1,ef[idc][6]),ec0);
                     transfer(f10);
                 }
                 warper->printEntity(f10,false);
                 // create new face
                 entity_t *f11 = warper->getFaceFromVertex(nvc,vs1,v1);
                 if (!f11)
                 {
                     f11 = warper->createFaceWithEdge(e1,t.get(1,ef[idc][7]),ec1);
                     transfer(f11);
                 }
                 warper->printEntity(f11,false);
                 // prepare storage
                 t1_t &face_related_t1_1 = hanging_data_t1.setData(*f1);
                 // find smallest address in between proc if ever for v0
                 const entity_t * v0s = smallestAddress(v0);
                 // find smallest address in between proc if ever for v1
                 const entity_t * v1s = smallestAddress(v1);

                 char frt1 = 2;
                 char frt2 = 1;
                 if (v0s < v1s) { frt1 = 1; frt2 = 2; }

                 face_related_t1_1[0] = warper->getEdgeFromVertex(nvc,vs1);
                 face_related_t1_1[frt1] = warper->getFaceFromVertex(nvc,vs1,v0);
                 face_related_t1_1[frt2] = warper->getFaceFromVertex(nvc,vs1,v1);

                 // Here we know that this face have to be splited for both proc if on proc bnd. We don't need to comunicate
                 // this event. We directely add this face to partition manager update treatement.
                 if (part_man->hasRemoteCopy(*f1))
                     edge_face_to_be_treated.insert(f1);

                 //====== Second face treatement
                 entity_t *f2 = t.get(2,ef[idc][1]);
                 assert(f2->size(0) == 3);

                 // collect information related to second face
                 collect(f2);

                 // create new edge
                 entity_t *e2 = warper->getEdgeFromVertex(nvc,vs2);
                 if (!e2)
                 {
                     e2 = warper->createEdge(nvc,vs2);
                     transfer(e2);
                 }
                 warper->printEntity(e2,false);
                 // create new face
                 entity_t *f20 = warper->getFaceFromVertex(nvc,vs2,v0);
                 if (!f20)
                 {
                     f20 = warper->createFaceWithEdge(e2,t.get(1,ef[idc][8]),ec0);
                     transfer(f20);
                 }
                 warper->printEntity(f20,false);
                 // create new face
                 entity_t *f21 = warper->getFaceFromVertex(nvc,vs2,v1);
                 if (!f21)
                 {
                     f21 = warper->createFaceWithEdge(e2,t.get(1,ef[idc][9]),ec1);
                     transfer(f21);
                 }
                 warper->printEntity(f21,false);
                 // prepare storage
                 t1_t &face_related_t1_t2 = hanging_data_t1.setData(*f2);

                 face_related_t1_t2[0] = warper->getEdgeFromVertex(nvc,vs2);
                 face_related_t1_t2[frt1] = warper->getFaceFromVertex(nvc,vs2,v0);
                 face_related_t1_t2[frt2] = warper->getFaceFromVertex(nvc,vs2,v1);

                 // Here we know that this face have to be splited for both proc if on proc bnd. We don't need to comunicate
                 // this event. We directely add this face to partition manager update treatement.
                 if (part_man->hasRemoteCopy(*f2))
                     edge_face_to_be_treated.insert(f2);

                 // collect information related to this tet
                 collect(&t);

                 //mid face creation
                 entity_t *f3 = warper->createFaceWithEdge(e1,t.get(1,ef[idc][10]),e2);
                 transfer(f3);

                 // tet creation
                 k = -1;
                 tet[++k] = warper->createTetWithFace(f11,f21,t.get(2,ef[idc][11]),f3);
                 //std::cout<<" 1nc "<<warper->tetAspectRatio(tet[k]);
                 tet[++k] = warper->createTetWithFace(f10,f3,t.get(2,ef[idc][12]),f20);
                 //std::cout<<" "<<warper->tetAspectRatio(tet[k]);
                 break;
             }
            // central node stategy
            default :
             {
                 newv = warper->createCOGVertex(t);
                 // uncoment to create a COG of all nodes used for tet (i.e. nodes of original tet + nodes of splited edges)
                 //for (i = 0; i < 4; i++) warper->accumulatePointCoord(t.get(0,i));
                 //newv = warper->createMeanPointCoord();

                 k = -1;
                 entity_t *ve[3];

                 // loop on faces
                 for (i = 0; i < 4; i++)
                 {
                     entity_t *f = t.get(2,i);
                     assert(f->size(0) == 3);
                     // collect information related to this face
                     collect(f);

                     // loop on face edges to set nodes
                     for (j = 0; j < 3; ++j)
                     {
                         v[j] = f->get(0,j);
                         entity_t **pv = hanging_down.getData(*f->get(1,j));
                         if (pv)
                         {
                             ve[j] = *pv;
                         }
                         else
                             ve[j] = nullptr;
                     }

                     // 3 edges split
                     if (ve[0] && ve[1] && ve[2])
                     {
                         // In this case and only in this one a face may be unsplit but all its edge are. In all other
                         // case its imposible to have a split face has not all edges are split.
                         // If face is not split we have to split it. By calling splitFace we do the test and split if needed.
                         // But we need to test here to know if partion manager have to be updated or not.
                         entity_t **pfc = hanging_down.getData(*f);
                         if (!pfc)
                         {
                             splitFace(*f);
                             // Here we know that this face have to be splited for both proc if on proc bnd. We don't need to comunicate
                             // this event. We directely add this face to partition manager update treatement.
                             if (part_man->hasRemoteCopy(*f))
                                 edge_face_to_be_treated.insert(f);
                         }

                         // create 4 tet from 4 sub face and newv
                         // no transfert (done by splitFace here or in splitTet)
                         tet[++k] = createTetFromVertex(newv,ve[0],ve[1],ve[2],trans);
                         tet[++k] = createTetFromVertex(newv,v[0],ve[0],ve[2],trans);
                         tet[++k] = createTetFromVertex(newv,v[1],ve[1],ve[0],trans);
                         tet[++k] = createTetFromVertex(newv,v[2],ve[2],ve[1],trans);
                     }
                     // first 2
                     else if (ve[0] && ve[1] )
                     {
                         // prepare storage
                         t2_t &face_related_t2 = hanging_data_t2.setData(*f);

                         // find smallest address in between proc if ever for ve[0]
                         const entity_t * ve0 = smallestAddress(ve[0]);
                         // find smallest address in between proc if ever for ve[1]
                         const entity_t * ve1 = smallestAddress(ve[1]);


                         // chose smallest address has split edge start
                         if (ve0 < ve1)
                         {
                             trans[Tve[1][2]] = trans[8] = true;
                             tet[++k] = createTetFromVertex(newv,ve[0],ve[1],v[1],trans);
                             tet[++k] = createTetFromVertex(newv,ve[0],v[2],v[0],trans);
                             tet[++k] = createTetFromVertex(newv,ve[0],v[2],ve[1],trans);
                             trans[Tve[1][2]] = trans[8] = false;
                             face_related_t2[0] = warper->getEdgeFromVertex(ve[0],v[2]);
                             face_related_t2[1] = warper->getEdgeFromVertex(ve[0],ve[1]);
                             face_related_t2[2] = warper->getFaceFromVertex(ve[0],ve[1],v[1]);
                             face_related_t2[3] = warper->getFaceFromVertex(ve[0],v[2],ve[1]);
                             face_related_t2[4] = warper->getFaceFromVertex(ve[0],v[2],v[0]);
                         }
                         else
                         {
                             trans[Tve[1][2]] = trans[8] = true;
                             tet[++k] = createTetFromVertex(newv,ve[1],ve[0],v[1],trans);
                             tet[++k] = createTetFromVertex(newv,ve[1],v[0],ve[0],trans);
                             tet[++k] = createTetFromVertex(newv,ve[1],v[0],v[2],trans);
                             trans[Tve[1][2]] = trans[8] = false;
                             face_related_t2[0] = warper->getEdgeFromVertex(ve[1],v[0]);
                             face_related_t2[1] = warper->getEdgeFromVertex(ve[1],ve[0]);
                             face_related_t2[2] = warper->getFaceFromVertex(ve[1],ve[0],v[1]);
                             face_related_t2[3] = warper->getFaceFromVertex(ve[1],v[0],ve[0]);
                             face_related_t2[4] = warper->getFaceFromVertex(ve[1],v[0],v[2]);
                         }

                         // Here we know that this face have to be splited for both proc if on proc bnd. We don't need to comunicate
                         // this event. We directely add this face to partition manager update treatement.
                         if (part_man->hasRemoteCopy(*f))
                             edge_face_to_be_treated.insert(f);
                     }
                     // second 2
                     else if (ve[2] && ve[1] )
                     {
                         // prepare storage
                         t2_t &face_related_t2 = hanging_data_t2.setData(*f);

                         // find smallest address in between proc if ever for ve[1]
                         const entity_t * ve1 = smallestAddress(ve[1]);
                         // find smallest address in between proc if ever for ve[2]
                         const entity_t * ve2 = smallestAddress(ve[2]);

                         // chose smallest address has split edge start
                         if (ve1 < ve2)
                         {
                             trans[Tve[1][2]] = trans[8] = true;
                             tet[++k] = createTetFromVertex(newv,ve[1],ve[2],v[2],trans);
                             tet[++k] = createTetFromVertex(newv,ve[1],v[0],v[1],trans);
                             tet[++k] = createTetFromVertex(newv,ve[1],v[0],ve[2],trans);
                             trans[Tve[1][2]] = trans[8] = false;
                             face_related_t2[0] = warper->getEdgeFromVertex(ve[1],v[0]);
                             face_related_t2[1] = warper->getEdgeFromVertex(ve[1],ve[2]);
                             face_related_t2[2] = warper->getFaceFromVertex(ve[1],ve[2],v[2]);
                             face_related_t2[3] = warper->getFaceFromVertex(ve[1],v[0],ve[2]);
                             face_related_t2[4] = warper->getFaceFromVertex(ve[1],v[0],v[1]);
                         }
                         else
                         {
                             trans[Tve[1][2]] = trans[8] = true;
                             tet[++k] = createTetFromVertex(newv,ve[2],ve[1],v[2],trans);
                             tet[++k] = createTetFromVertex(newv,ve[2],v[1],v[0],trans);
                             tet[++k] = createTetFromVertex(newv,ve[2],v[1],ve[1],trans);
                             trans[Tve[1][2]] = trans[8] = false;
                             face_related_t2[0] = warper->getEdgeFromVertex(ve[2],v[1]);
                             face_related_t2[1] = warper->getEdgeFromVertex(ve[2],ve[1]);
                             face_related_t2[2] = warper->getFaceFromVertex(ve[2],ve[1],v[2]);
                             face_related_t2[3] = warper->getFaceFromVertex(ve[2],v[1],ve[1]);
                             face_related_t2[4] = warper->getFaceFromVertex(ve[2],v[1],v[0]);
                         }
                         // Here we know that this face have to be splited for both proc if on proc bnd. We don't need to comunicate
                         // this event. We directely add this face to partition manager update treatement.
                         if (part_man->hasRemoteCopy(*f))
                             edge_face_to_be_treated.insert(f);
                     }
                     // third 2
                     else if (ve[2] && ve[0] )
                     {
                         // prepare storage
                         t2_t &face_related_t2 = hanging_data_t2.setData(*f);

                         // find smallest address in between proc if ever for ve[0]
                         const entity_t * ve0 = smallestAddress(ve[0]);
                         // find smallest address in between proc if ever for ve[2]
                         const entity_t * ve2 = smallestAddress(ve[2]);

                         // chose smallest address has split edge start
                         if (ve0 < ve2)
                         {
                             trans[Tve[1][2]] = trans[8] = true;
                             tet[++k] = createTetFromVertex(newv,ve[0],ve[2],v[0],trans);
                             tet[++k] = createTetFromVertex(newv,ve[0],v[2],v[1],trans);
                             tet[++k] = createTetFromVertex(newv,ve[0],v[2],ve[2],trans);
                             trans[Tve[1][2]] = trans[8] = false;
                             face_related_t2[0] = warper->getEdgeFromVertex(ve[0],v[2]);
                             face_related_t2[1] = warper->getEdgeFromVertex(ve[0],ve[2]);
                             face_related_t2[2] = warper->getFaceFromVertex(ve[0],ve[2],v[0]);
                             face_related_t2[3] = warper->getFaceFromVertex(ve[0],v[2],ve[2]);
                             face_related_t2[4] = warper->getFaceFromVertex(ve[0],v[2],v[1]);
                         }
                         else
                         {
                             trans[Tve[1][2]] = trans[8] = true;
                             tet[++k] = createTetFromVertex(newv,ve[2],ve[0],v[0],trans);
                             tet[++k] = createTetFromVertex(newv,ve[2],v[1],v[2],trans);
                             tet[++k] = createTetFromVertex(newv,ve[2],v[1],ve[0],trans);
                             trans[Tve[1][2]] = trans[8] = false;
                             face_related_t2[0] = warper->getEdgeFromVertex(ve[2],v[1]);
                             face_related_t2[1] = warper->getEdgeFromVertex(ve[2],ve[0]);
                             face_related_t2[2] = warper->getFaceFromVertex(ve[2],ve[0],v[0]);
                             face_related_t2[3] = warper->getFaceFromVertex(ve[2],v[1],ve[0]);
                             face_related_t2[4] = warper->getFaceFromVertex(ve[2],v[1],v[2]);
                         }
                         // Here we know that this face have to be splited for both proc if on proc bnd. We don't need to comunicate
                         // this event. We directely add this face to partition manager update treatement.
                         if (part_man->hasRemoteCopy(*f))
                             edge_face_to_be_treated.insert(f);
                     }
                     // first edge split
                     else if (ve[0] )
                     {
                         // prepare storage
                         t1_t &face_related_t1 = hanging_data_t1.setData(*f);

                         // edge ve[0] v[2] in the middle
                         trans[Tve[1][2]] = trans[8] = true;
                         tet[++k] = createTetFromVertex(newv,ve[0],v[2],v[0],trans);
                         //std::cout<<" 1c "<<warper->tetAspectRatio(tet[k]);
                         tet[++k] = createTetFromVertex(newv,ve[0],v[2],v[1],trans);
                         //std::cout<<" "<<warper->tetAspectRatio(tet[k]);
                         trans[Tve[1][2]] = trans[8] = false;
                         // find smallest address in between proc if ever for v[0]
                         const entity_t * v0 = smallestAddress(v[0]);
                         // find smallest address in between proc if ever for v[1]
                         const entity_t * v1 = smallestAddress(v[1]);

                         char frt1 = 2;
                         char frt2 = 1;
                         if (v0 < v1) { frt1 = 1; frt2 = 2; }

                         face_related_t1[0] = warper->getEdgeFromVertex(ve[0],v[2]);
                         face_related_t1[frt1] = warper->getFaceFromVertex(ve[0],v[2],v[0]);
                         face_related_t1[frt2] = warper->getFaceFromVertex(ve[0],v[2],v[1]);

                         // Here we know that this face have to be splited for both proc if on proc bnd. We don't need to comunicate
                         // this event. We directely add this face to partition manager update treatement.
                         if (part_man->hasRemoteCopy(*f))
                             edge_face_to_be_treated.insert(f);
                     }
                     // second edge split
                     else if (ve[1] )
                     {
                         // prepare storage
                         t1_t &face_related_t1 = hanging_data_t1.setData(*f);

                         // edge ve[1] v[0] in the middle
                         trans[Tve[1][2]] = trans[8] = true;
                         tet[++k] = createTetFromVertex(newv,ve[1],v[0],v[1],trans);
                         //std::cout<<" 2c "<<warper->tetAspectRatio(tet[k]);
                         tet[++k] = createTetFromVertex(newv,ve[1],v[0],v[2],trans);
                         //std::cout<<" "<<warper->tetAspectRatio(tet[k]);
                         trans[Tve[1][2]] = trans[8] = false;
                         // find smallest address in between proc if ever for v[1]
                         const entity_t * v1 = smallestAddress(v[1]);
                         // find smallest address in between proc if ever for v[2]
                         const entity_t * v2 = smallestAddress(v[2]);

                         char frt1 = 2;
                         char frt2 = 1;
                         if (v1 < v2) { frt1 = 1; frt2 = 2; }

                         face_related_t1[0] = warper->getEdgeFromVertex(ve[1],v[0]);
                         face_related_t1[frt1] = warper->getFaceFromVertex(ve[1],v[0],v[1]);
                         face_related_t1[frt2] = warper->getFaceFromVertex(ve[1],v[0],v[2]);

                         // Here we know that this face have to be splited for both proc if on proc bnd. We don't need to comunicate
                         // this event. We directely add this face to partition manager update treatement.
                         if (part_man->hasRemoteCopy(*f))
                             edge_face_to_be_treated.insert(f);
                     }
                     // third edge split
                     else if (ve[2] )
                     {
                         // prepare storage
                         t1_t &face_related_t1 = hanging_data_t1.setData(*f);

                         // edge ve[2] v[1] in the middle
                         trans[Tve[1][2]] = trans[8] = true;
                         tet[++k] = createTetFromVertex(newv,ve[2],v[1],v[0],trans);
                         //std::cout<<" 3c "<<warper->tetAspectRatio(tet[k]);
                         tet[++k] = createTetFromVertex(newv,ve[2],v[1],v[2],trans);
                         //std::cout<<" "<<warper->tetAspectRatio(tet[k]);
                         trans[Tve[1][2]] = trans[8] = false;
                         // find smallest address in between proc if ever for v[0]
                         const entity_t * v0 = smallestAddress(v[0]);
                         // find smallest address in between proc if ever for v[2]
                         const entity_t * v2 = smallestAddress(v[2]);

                         char frt1 = 2;
                         char frt2 = 1;
                         if (v0 < v2) { frt1 = 1; frt2 = 2; }

                         face_related_t1[0] = warper->getEdgeFromVertex(ve[2],v[1]);
                         face_related_t1[frt1] = warper->getFaceFromVertex(ve[2],v[1],v[0]);
                         face_related_t1[frt2] = warper->getFaceFromVertex(ve[2],v[1],v[2]);

                         // Here we know that this face have to be splited for both proc if on proc bnd. We don't need to comunicate
                         // this event. We directely add this face to partition manager update treatement.
                         if (part_man->hasRemoteCopy(*f))
                             edge_face_to_be_treated.insert(f);
                     }
                     // no cut
                     else
                     {
                         // no transfert
                         // no extra treatment
                         tet[++k] = createTetFromVertex(newv,v[0],v[1],v[2],trans);
                         //std::cout<<" nc "<<warper->tetAspectRatio(tet[k]);
                     }
                 }

                 // collect information related to this tet
                 collect(&t);
             }
        }
        //std::cout<<std::endl;

        // Transfer information to all newly created tet
        ++k;
        for (j = 0; j < k; ++j) transfer(tet[j]);

        if (newv)
        {
            // Transfer information from tet has  newv is inside this tet
            transfer(newv);

            // Transfer information from tet has every edge/face emanating from newv is inside this tet
            assert(newv->size(2));
            for (int i = 0,m = newv->size(2); i < m; ++i)
            {
                transfer(newv->get(2,i));
            }
            assert(newv->size(1));
            for (int i = 0,m = newv->size(1); i < m; ++i)
            {
                transfer(newv->get(1,i));
            }

        }

        // remove old tet
        warper->removeEntityWithAdj(t);

    }
    return newv;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
bool xMeshSplitUserCriteria < T,DM >::deleteFaceWithCondition(entity_t &face,entity_t &cface)
{
    // this face is related only to this tet : part boundary or hanging, proc boundary excluded
    if (face.size(target_dim) == 1 &&  !part_man->hasRemoteCopy(face))
    {

        // clean hanging information has this face has to be removed
        deleteHangingInfoForFace(face,cface);

        // remove face and its references
        warper->removeEntityWithAdj(face);

        return true;
    }
    else
        return false;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
void xMeshSplitUserCriteria < T,DM >::deleteHangingInfoForFace(entity_t &face,entity_t &cface)
{

    // clean hanging information for a face
    hanging_down.deleteData(face);
    hanging_up.deleteData(cface);
    // looping on edges of central face to found lateral faces
    int j,m;
    entity_t *fk;
    for (int k = 0; k < 3; ++k)
    {
        fk = nullptr;
        entity_t *edge = cface.get(1,k);
        assert(edge->size(2));
        for (j = 0,m = edge->size(2); j < m; ++j)
        {
            entity_t *lface = edge->get(2,j);
            if (lface != &cface)
            {
                entity_t **puface = hanging_up.getData(*lface);
                if (puface && *puface == &face)
                {
                    fk = lface;
                    break;
                }
            }
        }
        assert(fk);
        hanging_up.deleteData(*fk);
        hanging_up.deleteData(*edge);
    }

}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
bool xMeshSplitUserCriteria < T,DM >::deleteEdgeWithCondition(entity_t &edge,entity_t &midvertex)
{
    // this edge is related only to this entity (i.e. the one this edge belong to) :
    // part boundary or hanging, proc boundary excluded
    if (edge.size(target_dim) == 1 &&  !part_man->hasRemoteCopy(edge))
    {
        // clean hanging information has this edge has to be removed
        deleteHangingInfoForEdge(edge,midvertex);

        // remove edge and its references
        warper->removeEntityWithAdj(edge);

        return true;
    }
    else
        return false;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
void xMeshSplitUserCriteria < T,DM >::deleteHangingInfoForEdge(entity_t &edge,entity_t &midvertex)
{
    // loop on all surounding edge of central node to find
    // those who have hanging information to delete
    int l,n,k;
    for (l = 0,n = 0,k = midvertex.size(1); l < k && n < 2; ++l)
    {
        entity_t *se = midvertex.get(1,l);
        entity_t **h = hanging_up.getData(*se);
        if ( h != nullptr && *h == &edge )
        {
            hanging_up.deleteData(*se);
            ++n;
        }

    }
    assert(n == 2);
    hanging_down.deleteData(edge);
    hanging_up.deleteData(midvertex);
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
bool xMeshSplitUserCriteria < T,DM >::dontDeleteEntityWithCondition(entity_t &a,entity_t &b)
{
    return false;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
void xMeshSplitUserCriteria < T,DM >::hangingSynchro(entity_t &e)
{
    switch (e.getLevel())
    {
        case 1 :
         {
             splitEdge(e);
             break;
         }
        case 2 :
         {
             splitFace(e);
             break;
         }
        default :
         {
             throw -1456;
         }
    }
    edge_face_to_be_treated.insert(&e);
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
auto xMeshSplitUserCriteria < T,DM >::createTetFromVertex(entity_t *v0,entity_t *v1,entity_t *v2,entity_t *v3,bool *trans)->entity_t *
{
    entity_t *v[4];
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
    bool *transf = trans+6;

    // edges treatment
    entity_t *edge[6];
    for (int i = 0; i < 6; ++i)
    {
        entity_t *e = warper->getEdgeFromVertex(v[Tev[i][0]],v[Tev[i][1]]);
        if (!e)
        {
            e = warper->createEdge(v[Tev[i][0]],v[Tev[i][1]]);
            if (trans[i])
                transfer(e);
        }
        edge[i] = e;
    }

    // faces treatement
    entity_t *face[4];
    for (int i = 0; i < 4; ++i)
    {
        entity_t *f = warper->getFaceFromVertex(v[Tfv[i][0]],v[Tfv[i][1]],v[Tfv[i][2]]);
        if (!f)
        {
            f = warper->createFaceWithEdge(edge[Tfe[i][0]],edge[Tfe[i][1]],edge[Tfe[i][2]]);
            if (transf[i])
                transfer(f);
        }
        face[i] = f;
    }

    return warper->createTetWithFace(face[0],face[1],face[2],face[3]);
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
void xMeshSplitUserCriteria < T,DM >::removeDangling()
{
    // add all entity having  hanging  entity
    for (int i = 1; i < target_dim; ++i)
    {
        auto itle = warper->end(i);
        for (auto itl = warper->begin(i); itl != itle; ++itl)
        {
            entity_t *e = *itl;
            if (hanging_down.getData(*e) )
            {
                edge_face_to_be_treated.insert(e);
            }
            else if (i > 1)
            {
                if (hanging_data_t1.getData(*e))
                    edge_face_to_be_treated.insert(e);
                else if (hanging_data_t2.getData(*e))
                    edge_face_to_be_treated.insert(e);
            }

        }
    }

    // create container/manager
    xtool::xKeyContainerSendAndRecv < information_key_t > dangling_key_cont(*world);
    splitInfoManagerDanglingCount < entity_t, DM  > dangling_info;

    // set nb of connection : number of element related to entity
    for (auto e : edge_face_to_be_treated) dangling_info.setNbConnect(e,e->size(target_dim));

    // set comunication keys for bnd entity to be treated : first pass gathering on owner the maximum
    dangling_key_cont.accumulateKeysOwnerGather( edge_face_to_be_treated.begin(), edge_face_to_be_treated.end(), *kmsar );
    exchangeInformation(dangling_key_cont,dangling_info);

    // set comunication keys for bnd entity to be treated : second pass scattering maximum found to all
    dangling_key_cont.clearKeys();
    dangling_key_cont.accumulateKeysOwnerScatter( edge_face_to_be_treated.begin(), edge_face_to_be_treated.end(), *kmsar );
    exchangeInformation(dangling_key_cont,dangling_info);


    // treat dangling face based on max number of related entity
    typename splitInfoManagerDanglingCount < entity_t, DM  >::information_t nb;
    for (auto e : edge_face_to_be_treated)
    {
        // get max number of connected entity related to e
        dangling_info.getNbConnect(e,nb);
        // clear information what ever happend after. If below we delete e, data manager used in splitInfoManagerDanglingCount
        // will have problem. It have to be cleaned before removal.
        dangling_info.clearNbConnect(e);

        // if all duplicate entity on other proc are not connected to any element, or if no element is connect to local
        // entity, it is a dangling entity to be removed
        if (!nb)
        {

            entity_t **ppe = hanging_down.getData(*e);
            if (ppe)
            {
                entity_t *pe = *ppe;
                // face with hanging
                if (pe->getLevel() > 0)
                    // clean hanging information has this face has to be removed
                    deleteHangingInfoForFace(*e,*pe);
                // edge with hanging
                else
                    // clean hanging information has this edge has to be removed
                    deleteHangingInfoForEdge(*e,*pe);
            }
            else
            {
                // clean sub entity attached to e if any
                if (hanging_data_t1.getData(*e)) hanging_data_t1.deleteData(*e);
                else if (hanging_data_t2.getData(*e))
                    hanging_data_t2.deleteData(*e);
                else
                {
                    //wrong selection in edge_face_to_be_treated
                    throw -17832;
                }
            }

            // If on boundary remove from partiton manager
            // note : test is in clear method
            part_man->remove(*e);

            // remove face/edge and its references
            warper->removeEntityWithAdj(*e);
        }
    }

    return;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
void xMeshSplitUserCriteria < T,DM >::dontSelect(entity_t *e){}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
void xMeshSplitUserCriteria < T,DM >::doSelect(entity_t *e)
{
    if (part_man->hasRemoteCopy(*e))
        new_hanging_key_cont->accumulateKeys(*e,*kmso);
    return;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
void xMeshSplitUserCriteria < T,DM >::collect(entity_t *e)
{
    for (auto trans : trans_information_container) trans->collect(e);
    return;
}
//-------------------------------------------------------------------------------------------------
template < typename T, template < typename > class DM >
void xMeshSplitUserCriteria < T,DM >::transfer(entity_t *e)
{
    for (auto trans : trans_information_container) trans->transfer(e);
    return;
}
template < typename T, template < typename > class DM >
auto xMeshSplitUserCriteria < T,DM >::smallestAddress(entity_t *e)->const entity_t *
{
    // find smallest address in between proc if ever for given argument
    const entity_t * se = e;
    xtool::xConstPartitionObject < entity_t > po = part_man->getConstPartitionObject(*se);
    if (po.hasRemoteObject())
    {
        for (auto ro : po.getRemoteObjectsCollectionRange( ))
        {
            const entity_t * re = ro.getObjectAddress();
            if (se > re ) se = re;
        }
    }
    return se;
}

}                                                               // end namespace
#endif
