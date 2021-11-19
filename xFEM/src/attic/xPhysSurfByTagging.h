/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#ifndef ___XPHYSSURFBYTAGGING_H
#define ___XPHYSSURFBYTAGGING_H
#include <vector>

#include "xEntityToEntity.h"
#include "xRegion.h"
#include "xPhysSurfParameter.h"

// macro definition for tagging and appartenance methode
// note : macro defined as number for portability (no little/big indian nor storage size prb)
//
// for tag_entities
//
/// Entity status for "out"
//! <br/>note : STRICT and LOOS OUT status are related : LOOS_OUT_STAT=2*STRICT_OUT_STAT
//! don't change that as it is use somewhere
#define  STRICT_OUT_STAT      1 // 0000 0001
/// Entity status for "out touching iso-zero"
//! <br/>note : STRICT and LOOS OUT status are related : LOOS_OUT_STAT=2*STRICT_OUT_STAT
//! don't change that as it is use somewhere
#define  LOOS_OUT_STAT        2 // 0000 0010
/// Entity status for "included in iso-zero"
#define  ISO_ZERO_STAT        4 // 0000 0100
/// Entity status for "cut by iso-zero"
#define  CUT_STAT             8 // 0000 1000
/// Entity status for "out"
//! <br/>note : STRICT and LOOS IN status are related : LOOS_IN_STAT=2*STRICT_IN_STAT
//! don't change that as it is use somewhere
#define  STRICT_IN_STAT      16 // 0001 0000
/// Entity status for "in touching iso-zero"
//! <br/>note : STRICT and LOOS IN status are related : LOOS_IN_STAT=2*STRICT_IN_STAT
//! don't change that as it is use somewhere
#define  LOOS_IN_STAT        32 // 0010 0000

// for tag_support
//
/// support status for "cut by iso-zero Elt wise"
#define  SCUTEW_STAT          1 // 0000 0001
/// support status for "cut by iso-zero"
#define  SCUT_STAT            2 // 0000 0010
/// support status for "in and touched by iso-zero"
#define  SREL_IN_STAT         4 // 0000 0100
/// support status for "out and touched by iso-zero"
#define  SREL_OUT_STAT        8 // 0000 1000

// for both
//
/// no status for entity and support tag
#define   NO_STAT             0 // 0000 0000

// mask for internal purpose
#define TS_RELATED           12 // 0000 1100 =  SREL_IN_STAT || SREL_OUT_STAT
#define TE_LOOS              34 // 0010 0010 =  LOOS_OUT_STAT || LOOS_IN_STAT

// mask for appartenance methode
// 
/// mask for supportCoversIn
//! <br/> for tag_entities. It is equivalent to test LOOS_IN_STAT || ISO_ZERO_STAT || CUT_STAT || STRICT_IN_STAT
#define TE_INCOVER          60 // 0011 1100 
//
/// mask for supportCoversOut
//! <br/> for tag_entities. It is equivalent to test LOOS_OUT_STAT || ISO_ZERO_STAT || CUT_STAT || STRICT_OUT_STAT
#define TE_OUTCOVER         15 // 0000 1111 
//
/// mask for supportCutStrictly
//! <br/> for tag_support. It is equivalent to test SCUTEW_STAT || SCUT_STAT
#define TS_CUTSTRICTLY       3 // 0000 0011
//
/// mask for supportBoundary
//! <br/> for tag_support. It is equivalent to test SCUTEW_STAT || SCUT_STAT || SREL_IN_STAT || SREL_OUT_STAT
#define TS_BOUNDARY         15 // 0000 1111 
//
/// mask for strictIn
//! <br/> for tag_entities. It is equivalent to test LOOS_IN_STAT || STRICT_IN_STAT
#define TE_INDOMAIN         48 // 0011 0000 
//
/// mask for coversIn
//! <br/> for tag_entities. It is equivalent to test LOOS_IN_STAT || STRICT_IN_STAT || CUT_STAT
#define TE_INDOMAINCUT      56 // 0011 1000 
//
/// mask for strictOut
//! <br/> for tag_entities. It is equivalent to test LOOS_OUT_STAT || STRICT_OUT_STAT
#define TE_OUTDOMAIN          3 // 0000 0011
//
/// mask for coversOut
//! <br/> for tag_entities. It is equivalent to test LOOS_OUT_STAT || STRICT_OUT_STAT || CUT_STAT
#define TE_OUTDOMAINCUT      11 // 0000 1011 


// uncoment to get time statistique
//#define TIMING_MONITORING 1
#ifdef TIMING_MONITORING
#include "xChrono.h"
#endif

namespace xfem
{

class xMesh;
class xLevelSet;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xcut::xPhysSurfByTagging class ///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// html for doxygen : please respect html tag will changing  comments
/// This class cut, tag and classify a mesh according to the position of the iso-zero of a level-set.
//! <br/> 
//! <br/>         PREREQUIRED : 
//! <ul><li>      in 2D, mesh must have consistent element orientation (same normal between two adjacent elements, in particular when xcut::xPhysSurfByTagging is recursive) <li/>
//! <ul/>
//! <br/> 
//! <br/>         Let put in place some terminological definition :
//! <br/>         Given the entity e, the support is defined as the set of elements conected to the entity e.
//! <br/>         Let Nse be the set of nodes of the support of "e".
//! <br/>         Let Ne  be the set of nodes of the entity "e".
//! <br/>         The folowing definition are used :
//! <table border="1" cellpadding="2" cellspacing="2" width="70%" style="text-align:center">
//! <tr ><td>     Math         </td><td> Definition </td></tr>
//! <tr><td>     smin          </td><td>   It is the minimun value of the level-set compute on nodes of Nse       </td></tr>
//! <tr><td>     smax          </td><td>   It is the maximum value of the level-set compute on nodes of Nse       </td></tr>
//! <tr><td>     min           </td><td>   It is the minimun value of the level-set compute on nodes of Ne        </td></tr>
//! <tr><td>     max           </td><td>   It is the maximun value of the level-set compute on nodes of Ne        </td></tr>
//
//!</table>
//! <br/> 
//! <br/>         In this class many methodes and concept use strings "in" and "out" to refer repectively to negative and positive value of the level set.
//! <br/> 
//! <br/>         All the work of cutting, tagging and classify is done at construction time.
//! <br/>         The only methodes wich may modify this inintial work is update. All others are "read only" methodes.
//! <br/> 
//! <br/>         The number of xcut::xPhysSurfByTagging used  at the same time is not limited but the way their iso-zero interact is.
//!               The rule of tombe is the folowing "a element of the mesh can't be cut more then 2 times"
//! <br/> 
//! <br/>         The user have a large choice of methodes to check placement of entities or their support compared to position
//!               of the iso-zero of the level-set.
//! <br/>         The folowing table gives all those methodes and their mathematical signification  :
//! <br/> 
//  ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//! <table border="1" cellpadding="2" cellspacing="2" width="100%" style="text-align:center">
//! <tr ><td>     Name                     </td><td> Atomic </td><td>       Math        </td><td>                        included in                                  </td></tr>
//! <tr ><td>  supportCoversIn             </td><td>        </td><td>     smin<0        </td><td>                                                                     </td></tr>
//! <tr ><td>  supportCoversOut            </td><td>        </td><td>     smax>0        </td><td>                                                                     </td></tr>
//! <tr ><td>  supportBoundary             </td><td>        </td><td>   smin.smax<=0    </td><td>                                                                     </td></tr>
//! <tr ><td> supportCutStrictly           </td><td>        </td><td>    smin.smax<0    </td><td>     supportBoundary,supportCoversIn,supportCoversout                </td></tr>
//! <tr ><td> supportCutStrictlyEltWise    </td><td>   X    </td><td>  smin.smax<0 & C  </td><td> supportCutStrictly,supportBoundary,supportCoversIn,supportCoversout </td></tr>
//! <tr ><td> noTouchingIn                 </td><td>   X    </td><td>     max<0         </td><td>                 strictIn,coversIn                                   </td></tr>
//! <tr ><td> touchingIn                   </td><td>   X    </td><td>  min<0 & max=0    </td><td>                 strictIn,coversIn                                   </td></tr>
//! <tr ><td> inIsoZero                    </td><td>   X    </td><td>    min=max=0      </td><td>                                                                     </td></tr>
//! <tr ><td> cutStrictly                  </td><td>   X    </td><td>    min.max<0      </td><td>                 coversIn,coversOut                                  </td></tr>
//! <tr ><td> touchingOut                  </td><td>   X    </td><td>  min=0 & max>0    </td><td>                 strictOut,coversOut                                 </td></tr>
//! <tr ><td> noTouchingOut                </td><td>   X    </td><td>     min>0         </td><td>                 strictOut,coversOut                                 </td></tr>
//! <tr ><td> strictIn                     </td><td>        </td><td>  min<0 & max<=0   </td><td>                                                                     </td></tr>
//! <tr ><td> coversIn                     </td><td>        </td><td>     min<0         </td><td>                                                                     </td></tr>
//! <tr ><td> strictOut                    </td><td>        </td><td>  min>=0 & max>0   </td><td>                                                                     </td></tr>
//! <tr ><td> coversOut                    </td><td>        </td><td>     max>0         </td><td>                                                                     </td></tr>
//! </table>
//  ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//
//! <br/>         where :
//! <br/>             - C is a additional condition wich is true if at least one of the element of the support is stricly cut by the iso-zero of the level set (min.max<0)
//! <br/>             - a X in Atomic column means that this methodes reflect directely the status of the  entity or support tag attached to the entity tested
//! <br/>             - column "included in" gives methodes that include appartence methode of the ligne considered. Including here, means containing at least the set of entity 
//!                     that will respond true to appartence methode of the ligne considered.
//! <br/>         
//! <br/>        The methode dealing with support have only a meanning for entity of the mesh. On the over hand all other appartenance methode will give
//!              a information for entity of the mesh and entity of the submeshs generated by this class.
//! <br/>         
//! <br/>        When using xIntegrationRulePartition remember that it is the entity of the deepest level of recurtion  wich is filtered. Using cutStrictly as filter will then give
//!              no treatement (output, integration, ....) has a sumesh entity is not cut by iso-zero. By constuction it may only be touching iso-zero or strictly in a domain.
//!              The same filter used with xIntegrationRuleBasic will give all the element of the mesh cutted strictly by the iso-zero.
//
class xcut::xPhysSurfByTagging
{
    public:
        //
        /// The constructor is doing most of the job of this class
        xcut::xPhysSurfByTagging(xLevelSet & ls_, xcut::xPhysSurfParameter param = xcut::xPhysSurfParameter() );
        virtual ~xcut::xPhysSurfByTagging();

        // public methodes ///////////////

        /// Get the mesh 
        xMesh &getMesh();
        /// Get the mesh (const version)
        const xMesh &getMesh() const;

        /// Get a pointeur to the mesh 
        xMesh * getMeshPtr();
        /// Get a pointeur to the mesh (const version)
        const xMesh * getMeshPtr() const;

        /// Get the iso-zero mesh
        xMesh *  getMesh_bnd();
        /// Get the iso-zero mesh (const version)
        const xMesh *  getMesh_bnd() const;

        //
        /// Support appartenance methode : return true if entity "e" verify smin < 0
        bool supportCoversIn(AOMD::mEntity* e) const;

        //
        /// Support appartenance methode : return true if entity "e" verify smax > 0
        bool supportCoversOut(AOMD::mEntity* e) const;

        /// Support appartenance methode : return true if entity "e" verify smin.smax <= 0
        bool supportBoundary(AOMD::mEntity* e) const; 
        //
        /// Support appartenance methode : return true if entity "e" verify smin.smax < 0
        bool supportCutStrictly(AOMD::mEntity* e) const;
        //
        /// Support appartenance methode : return true if entity "e" verify smin.smax < 0 and at least one of the element of the support verify min.max<0
        bool supportCutStrictlyEltWise(AOMD::mEntity* e) const;
        
        /// Entity atomic appartenance methode : return true if entity "e" verify max < 0
        bool noTouchingIn(AOMD::mEntity* e) const;
        //
        /// Entity atomic appartenance methode : return true if entity  "e" verify  min<0 & max=0
        bool touchingIn(AOMD::mEntity* e) const;
        //
        /// Entity atomic appartenance methode : return true if entity  "e" verify  min=max=0
        bool inIsoZero(AOMD::mEntity* e) const;
        //
        /// Entity atomic appartenance methode : return true if entity  "e" verify  min.max<0
        bool cutStrictly(AOMD::mEntity* e) const;
        //
        /// Entity atomic appartenance methode : return true if entity  "e" verify min=0 & max>0
        bool touchingOut(AOMD::mEntity* e) const;
        //
        /// Entity atomic appartenance methode : return true if entity  "e" verify  min>0
        bool noTouchingOut(AOMD::mEntity* e) const;

        //
        /// Entity appartenance methode : return true if entity  "e" verify  min<0 & max<=0
        bool strictIn(AOMD::mEntity* e) const;
        //
        /// Appartenance methode : return true if entity  "e" verify  min<0
        bool coversIn(AOMD::mEntity* e) const;

        //
        /// Appartenance methode : return true if entity  "e" verify  min>=0 & max>0
        bool strictOut(AOMD::mEntity* e) const;
        //
        /// Appartenance methode : return true if entity  "e" verify  max>0
        bool coversOut(AOMD::mEntity* e) const;

        /// Give dimension of the mesh attached to this xcut::xPhysSurfByTagging
        int  dim() const;

        /// Give the level set attached to this xcut::xPhysSurfByTagging
        xLevelSet & getLevelSet() { return ls; }
        /// Give the level set attached to this xcut::xPhysSurfByTagging (const version)
        const xLevelSet & getLevelSet() const { return ls; }

        /// Give classifyer ("in") attached to this xcut::xPhysSurfByTagging
        xEntityToEntity &getClassifyerIn(){return classify_in; }
        /// Give classifyer ("out") attached to this xcut::xPhysSurfByTagging 
        xEntityToEntity &getClassifyerOut(){return classify_out; }

        /// Update classification and cutting of the domaine (when level set attached to this xcut::xPhysSurfByTagging change)
        void  update(bool fit = true, bool keep_old_partition = false,  bool recursive_ = false);

        /// for promotion xcut::xRefCutToAOMD have to be a friend of xcut::xPhysSurfByTagging
        friend class xcut::xRefCutToAOMD;
	friend class xPhysSurfVLS;
	friend class xVLSTriangleCutAttachableData;
	

#ifdef PARALLEL
        // In construct_ for paralle version a "exchange" of tag value  have to be done at the end to complete tagging

        // nota : set public because they are called as calback of DataExchange and have then to be accesible from their

        // Use of DataExchange function to "exchange"  tag value
        // Prerequist of DataExchange:
        typedef short data_type; // by default the 2 tag values are exchange for a entity betwen process. Here a short is considered to be 2 bytes long (2 char of 1 bytes long)
        // The send methode fill it's second argument with tag value from entity of it's first argument
        void send(AOMD::mEntity *, data_type &) const;
        // The receive methode deduce from it's second argument wich tag value it have to set (in conjonction of existing ones) for entity of it's first argument
        void receive(AOMD::mEntity *, const std::vector < data_type > &) const;

#endif
        // public members ////////////////

    protected:
        // private members ////////////////
        /// pointer to mesh given by level-set
        xMesh*  mesh;
        /// pointer to mesh of the iso-zero of the level-set constructed by this tag
        xMesh*  mesh_bnd;
        /// region constructed from mesh given by level-set
        xRegion region;

        /// classifyer for "in domain" according to the sens given by old xcut::xPhysSurf
        xEntityToEntity classify_in;
        /// classifyer for "out domain" according to the sens given by old xcut::xPhysSurf
        xEntityToEntity classify_out;

        /// level set reference
        xLevelSet & ls;

        /// boolen to do or not the fitting of the level-set according to fittol values
        bool fit;
        
        /// tolerance for the fitting of the level-set
        double fittol;

        /// Only create higher dimension partitions (only sub tets in 3D, sub_tri in 2D)
        bool only_higher_dim_partition;

        // MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---
        // for now this are considered as semi public tag but for MULTILEVELSET they must be considered
        // as temporary
        // MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---MULTILEVELSET---
 
        // html for doxygen : please respect html tag will changing this comment

        ///               tags needed for appartenance methode response
        //!        
        //!        
        //! <br/>         This tag of the entity (or sub entity of the cutting process) is designed to define the status of
        //! <br/>         a entitie against its position from iso-zero surface of the level-set.
        //! <br/>         The folowing table decribe those status, theire matematical meaning and entity capabilty against such status.
        //
        //! <table border="1" cellpadding="2" cellspacing="2" width="100%" style="text-align:center">
        //
        //  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr ><td> Entity statut           </td><td>    Math status    </td><td> status value macro </td><td> Vertex </td><td>   Edge  </td><td>  Face </td><td> Volume </td></tr>
        //  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr><td>       out                </td><td>      min>0        </td><td>  STRICT_OUT_STAT   </td><td>   X    </td><td>    X    </td><td>   X   </td><td>   X    </td></tr>
        //  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr><td>  out touching iso-zero   </td><td>  min>=0 & max>0   </td><td>   LOOS_OUT_STAT    </td><td>        </td><td>    X    </td><td>   X   </td><td>   X    </td></tr>
        //  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr><td>  included in iso-zero    </td><td>     min=max=0     </td><td>   ISO_ZERO_STAT    </td><td>   X    </td><td> X(2D,3D)</td><td> X(3D) </td><td>        </td></tr>
        //  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr><td>    cut by iso-zero       </td><td>  min<0 && max>0   </td><td>       CUT_STAT      </td><td>        </td><td>    X   </td><td>   X   </td><td>   X    </td></tr>
        //  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr><td>  in touching iso-zero    </td><td>  max<=0 & min<0   </td><td>   LOOS_IN_STAT     </td><td>        </td><td>    X    </td><td>   X   </td><td>   X    </td></tr>
        //  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr><td>       in                 </td><td>      max<0        </td><td>  STRICT_IN_STAT    </td><td>   X    </td><td>    X    </td><td>   X   </td><td>   X    </td></tr>
        //  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //
        //!</table>
        //
        //! <br/>         where :
        //! <br/>               in capabilty column X indicated that this status is possible for this kind of entity
        //! <br/>               in capabilty column void indicated that this status is impossible for this kind of entity
        //! <br/>         
        //! <br/>         All these status are mutualy exclusive (in a capabilty column only one possible ligne may be set)
        //! <br/>         Nota : a over status not mentioned in the above table is relateded to unset tag : NO_STAT 
        unsigned int tag_entities;

        ///               tags needed for appartenance methode response
        //!        
        //!        
        //! <br/>         This tag of the support of a entity  is designed to define the status of
        //! <br/>         the support of a entity against it's position from iso-zero surface of the level-set.
        //! <br/>         The folowing table decribe those status, theire matematical meaning and entity capabilty against such status.
        //
        //! <table border="1" cellpadding="2" cellspacing="2" width="100%" style="text-align:center">
        //
        //  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr><td>   Support statut            </td><td>    Math status      </td><td> status value macro </td><td> Vertex  </td><td> Edge </td><td> Face </td><td> Volume </td></tr>
        //  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr><td> in and touched by iso-zero  </td><td>      smax = 0       </td><td>    SREL_IN_STAT    </td><td>   X    </td><td>   X   </td><td>  X  </td><td>   X    </td></tr>
        //  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr><td> out and touched by iso-zero </td><td>      smin = 0       </td><td>    SREL_OUT_STAT   </td><td>   X    </td><td>   X   </td><td>  X  </td><td>   X    </td></tr>
        //  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr><td>    cut by iso-zero          </td><td>  smin.smax <0 && !C </td><td>      SCUT_STAT     </td><td>   X    </td><td>   X   </td><td> X(3D)</td><td>        </td></tr>
        //  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr><td> cut by iso-zero Elt wise    </td><td>  smin.smax <0 && C  </td><td>      SCUTEW_STAT   </td><td>   X    </td><td>   X   </td><td>  X   </td><td>   X    </td></tr>
        //  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //! <tr><td>      not cut                </td><td>  smax<0  || smin>0  </td><td>        NO_STAT     </td><td>   X    </td><td>   X   </td><td>  X   </td><td>   X    </td></tr>
        //  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //
        //!</table>
        //
        //! <br/>         where :
        //! <br/>               C is a additional condition wich is true if at least one of the element of the support is stricly cut by the iso-zero of the level set (min.max<0)
        //! <br/>               in capabilty column X indicated that this status is possible for this kind of entity
        //! <br/>               in capabilty column void indicated that this status is impossible for this kind of entity
        //! <br/>         
        //! <br/>         All these status are mutualy exclusive (in a capabilty column only one possible ligne may be set).
        unsigned int tag_support;


        /// code for appartenance methode
        //! thise code are compare in appartenance methode to respond true or false
        //! they are constructed from macro (defined above) in the constructor for portability reason (big/little indian ...)
        unsigned short int code_support_cover_in;
        unsigned short int code_support_cover_out;

        /// Promotors container
        //! this container give the ability to tag new sub entity with other xcut::xPhysSurf/level-set tags
        //!<br/> it store the other xcut::xPhysSurf to be able to acces to theire promotion methode 
        std::vector<xcut::xPhysSurfByTagging *> promotors; 
        /// Promotor status for this level set
        char promotor_status; 

        // private methodes ///////////////

        /// create tagging , cut and classify
        void construct_(bool fit = true, bool keep_old_partition_flag = false, bool recursive = false);
 
        /// Promotion methode to set Promotor status according to the tag of the entity given
        void  setPromotorStatus(AOMD::mEntity *);

        /// Promotion :
        //! Promotion itself Input entity will receive the same tag as the last entity passed to setPromotorStatus.
        void  promoteStatus(AOMD::mEntity *);
	/// tag eeventually attached by a previous call to promoteStatus is removed.
	void  unPromoteStatus(AOMD::mEntity *);

};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xcut::xPhysSurfByTagging class ///////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// xcut::xPhysSurfByTaggingException class //////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// interface derived class of standart exception for xcut::xPhysSurfByTagging
class xcut::xPhysSurfByTaggingException : public std::exception
{
    public:
        xcut::xPhysSurfByTaggingException(std::string,std::string,int,std::string,std::string);
        ~xcut::xPhysSurfByTaggingException() throw( );
        virtual const char * what() const throw( );
        // public members ////////////////

    private:
        // private members ////////////////

        std::string msg;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// End xcut::xPhysSurfByTaggingException  class /////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

} // end of namespace

#endif
