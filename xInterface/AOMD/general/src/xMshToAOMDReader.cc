#include <unordered_map>
#include <functional>
#include <cstring>
#include "xMshToAOMDReader.h"
#include "xReaderException.h"
#include "xAttachedDataManagerAOMD.h"
#include "mMesh.h"
#include "mEntity.h"
#include "mTet.h"
#include "mHex.h"
#include "mPrism.h"
#include "mFace.h"
#include "mEdge.h"
#include "mVertex.h"

using namespace std;
using namespace AOMD;
using namespace xinterface::aomd;

namespace xmeshtool
{
mEntity * genELM(mMesh &mesh, mVertex *nod[],int iTyp,int iGrp)
{
    mEntity *theEntity = nullptr;

    switch (iTyp)
    {
        case 4 :
            theEntity = static_cast < mEntity * >( mesh.createTetWithVertices(nod[0],nod[1],nod[2],nod[3],mesh.getGEntity(iGrp,3)));
            break;
        case 5 :
            theEntity = static_cast < mEntity * >( mesh.createHexWithVertices(nod[0],nod[1],nod[2],nod[3],nod[4],nod[5],nod[6],nod[7],
                                                                              mesh.getGEntity(iGrp,3)));
            break;
        case 6 :
            theEntity = static_cast < mEntity * >( mesh.createPrismWithVertices(nod[0],nod[1],nod[2],nod[3],nod[4],nod[5],
                                                                                mesh.getGEntity(iGrp,3)));
            break;
        case 2 :
            theEntity = static_cast < mEntity * >( mesh.createFaceWithVertices(nod[0],nod[1],nod[2],mesh.getGEntity(iGrp,2)));
            break;
        case 3 :
            theEntity = static_cast < mEntity * >( mesh.createFaceWithVertices(nod[0],nod[1],nod[2],nod[3],mesh.getGEntity(iGrp,2)));
            break;
        case 1 :
         {
             theEntity = static_cast < mEntity * >( mesh.createEdge(nod[0],nod[1],mesh.getGEntity(iGrp,1)));
             break;
         }
        case 15 :
         {
             mVertex *v = nod[0];
             {
                 v->classify(mesh.getGEntity(iGrp,0));
             }
             theEntity = static_cast < mEntity * >( v );
             break;
         }
        default :
         {
             std::ostringstream oss;
             oss << "Reader do not know this type of element: "<<iTyp<<std::endl;
             throw xReaderException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
         }
    }
    return theEntity;
}
void readNode(istream &f,mMesh &mesh, mVertex * nod[],int iElm,int iNbNod)
{
    int id;
    for (int i = 0; i < iNbNod; i++)
    {
        f >> id;
        nod[i] = mesh.getVertex(id);
        if (!nod[i])
        {
            std::ostringstream oss;
            oss << " Unknown vertex id "<<id<<" while scaning "<<i<<"th node of element "<<iElm<<std::endl;
            throw xReaderException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
        }
    }
}

void endBlock(istream &f,const char *s)
{
    char line[256];
    f.getline (line,256); // end line
    f.getline (line,256); // new line
    if (f.eof() || strncmp(&line[0],s,strlen(s)))
    {
        std::ostringstream oss;
        oss << s<<" is missing"<<std::endl;
        throw xReaderException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
    }
    return;
}
void readBlockMeshFormat(istream &f,double & gmshFormat)
{
    int i1,i2;
    f >>  gmshFormat >> i1 >> i2;
    cout << "gmsh mesh format " << gmshFormat << endl;
    endBlock(f,"$EndMeshFormat");
    return;
}

void readBlockNode(istream &f,mMesh &mesh,const char *s)
{
    int NbNod;
    f >> NbNod;
    for (int i = 0; i < NbNod; i++)
    {
        int iNod;
        double x,y,z;
        f >> iNod >> x >> y >> z;
        mesh.createVertex(iNod,x,y,z,nullptr);
    }
    endBlock(f,s);
}
template < bool read_id >
void elemCreate(mMesh &mesh, mVertex *nod[],int iElm, int iTyp,int iGrp,xAttachedDataManagerAOMD < int > &entities_id);
template < >
void elemCreate < true >(mMesh &mesh, mVertex *nod[],int iElm, int iTyp,int iGrp,xAttachedDataManagerAOMD < int > & entities_id)
{
    mEntity* elem = genELM(mesh, nod,iTyp, iGrp);
    entities_id.setData(*elem) = iElm;
}
template < >
void elemCreate < false >(mMesh &mesh, mVertex *nod[],int iElm, int iTyp,int iGrp,xAttachedDataManagerAOMD < int > & entities_id)
{
    genELM(mesh, nod,iTyp, iGrp);
}

template < bool read_id >
void readBlockELM(istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > & entities_id)
{
    int NbElm,iElm, iNbNod,iTyp,iGrp,iNbSub;
    mVertex *nod[100];
    f >> NbElm;
    for (int i = 0; i < NbElm; ++i)
    {
        // element type identification
        f >> iElm >> iTyp >> iGrp >> iNbSub >> iNbNod;

        // reading node
        readNode(f,mesh,nod,iElm,iNbNod);

        // element creation
        elemCreate < read_id >(mesh,nod,iElm,iTyp,iGrp,entities_id);
    }
    endBlock(f,"$ENDELM");
}
template < bool read_id >
void readBlockElements(istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > & entities_id)
{
    int NbElm,iElm,iParDom, iNbNod,iTyp,iGrp,iNbSub;
    mVertex *nod[100];
    f >> NbElm;
    for (int i = 0; i < NbElm; ++i )
    {
        f >> iElm >> iTyp >> iParDom >> iGrp >> iNbSub;

        switch (iTyp)
        {
            case 4 :
                iNbNod = 4;
                break;
            case 5 :
                iNbNod = 8;
                break;
            case 6 :
                iNbNod = 6;
                break;
            case 2 :
                iNbNod = 3;
                break;
            case 3 :
                iNbNod = 4;
                break;
            case 1 :
                iNbNod = 2;
                break;
            case 15 :
                iNbNod = 1;
                break;
            default :
             {
                 ostringstream oss;
                 oss << "Reader do not know this type of element: "<<iTyp<<endl;
                 throw xReaderException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
             }
        }

        // reading node
        readNode(f,mesh,nod,iElm,iNbNod);

        // element creation
        elemCreate < read_id >(mesh,nod,iElm,iTyp,iGrp,entities_id);
    }
    endBlock(f,"$EndElements");

}
typedef std::unordered_map < std::string,std::function < void ( istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > & entities_id ) > > tasks_t;

void generalLoop(istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > & entities_id,tasks_t &tasks)
{
    char line[256];
    while (1)
    {

        // seek block
        do
        {
            f.getline (line,256);
            if (f.eof()) break;
        } while (line[0] != '$');

        // no more block => end
        if (f.eof()) break;

        // search for coded block
        auto itf = tasks.find(line);

        //if unknown block
        if (itf == tasks.end())
        {
            // read unknow block till next token or eof
            cout<<"Block "<<line<<" skiped"<<endl;
            do
            {
                f.getline (line,256);
                if (f.eof()) break;
            } while (line[0] != '$');
            if (f.eof()) break;     // stop infinite loop
        }
        else
            itf->second(f,mesh,entities_id);  // treat known block


    }

    // classify unclassified vertices  ??? what for ???
    for (int i = 1; i < 4; i++)
    {
        if (mesh.size(i))
        {
            for (mMesh::iterall it = mesh.beginall(i); it != mesh.endall(i); ++it)
            {
                mEntity *e = *it;
                if (!e->getClassification())
                    e->print();
                for (int j = 0; j < e->size(0); j++)
                {
                    mEntity *v = e->get(0,j);
                    //	      v->print();
                    if (!v->getClassification() && e->getClassification())
                        v->classify(e->getClassification());
                }
            }
        }
    }

}

template < >
void xMshReader(istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > & entities_id)
{
    double gmshFormat = 0.;
    tasks_t tasks;
    tasks.insert(std::make_pair("$MeshFormat",[&gmshFormat](istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > &entities_id) {
                                    readBlockMeshFormat(f,gmshFormat);
                                }));
    tasks.insert(std::make_pair("$Nodes",[] ( istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > &entities_id ) {
                                    readBlockNode(f,mesh,"$EndNodes");
                                }));
    tasks.insert(std::make_pair("$NOD",[] ( istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > &entities_id ) {
                                    readBlockNode(f,mesh,"$ENDNOD");
                                }));
    tasks.insert(std::make_pair("$ELM",[] ( istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > &entities_id ) {
                                    readBlockELM < true >(f,mesh,entities_id);
                                }));
    tasks.insert(std::make_pair("$Elements",[&gmshFormat](istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > &entities_id) {
                                    if (gmshFormat > 2)
                                        readBlockElements < true >(f,mesh,entities_id);
                                    else
                                    {
                                        std::ostringstream oss;
                                        oss <<"Mesh format is not 2 or above will reading 'Elements' block available only in those versions!"<<std::endl;
                                        throw xReaderException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
                                    }
                                }));
    generalLoop(f,mesh,entities_id,tasks);

}

template < >
void xMshReader(istream &f,mMesh &mesh)
{
    xAttachedDataManagerAOMD < int > entities_id;
    double gmshFormat = 0.;
    tasks_t tasks;
    tasks.insert(std::make_pair("$MeshFormat",[&gmshFormat](istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > &entities_id) {
                                    readBlockMeshFormat(f,gmshFormat);
                                }));
    tasks.insert(std::make_pair("$Nodes",[] ( istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > &entities_id ) {
                                    readBlockNode(f,mesh,"$EndNodes");
                                }));
    tasks.insert(std::make_pair("$NOD",[] ( istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > &entities_id ) {
                                    readBlockNode(f,mesh,"$ENDNOD");
                                }));
    tasks.insert(std::make_pair("$ELM",[] ( istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > &entities_id ) {
                                    readBlockELM < false >(f,mesh,entities_id);
                                }));
    tasks.insert(std::make_pair("$Elements",[&gmshFormat](istream &f,mMesh &mesh, xAttachedDataManagerAOMD < int > &entities_id) {
                                    if (gmshFormat > 2)
                                        readBlockElements < false >(f,mesh,entities_id);
                                    else
                                    {
                                        std::ostringstream oss;
                                        oss <<"Mesh format is not 2 or above will reading 'Elements' block available only in those versions!"<<std::endl;
                                        throw xReaderException(oss.str(),__FILE__,__LINE__,__DATE__,__TIME__);
                                    }
                                }));
    generalLoop(f,mesh,entities_id,tasks);

}
}            // end namespace

