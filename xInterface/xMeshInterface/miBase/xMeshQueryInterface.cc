#include <xEntity.h>
#include <xEntityIterators.h>
#include <xMeshQueryInterface.h>
#include <xStringManager.h>

namespace xinterface {

namespace xmeshinterface {

    using std::ofstream;
    using std::string;

    xMeshQueryInterface::xMeshQueryInterface() {};

    unsigned int xMeshQueryInterface::lookupQueryDataId( const string& s )
    {
        return QueryStringManager.getId( s );
    };

    void xMeshQueryInterface::setName( const string _name )
    {
        name = _name;
        tag = QueryStringManager.getId( _name );
    };

    string xMeshQueryInterface::getName() const { return name; };

    xEntity xMeshQueryInterface::getGhostXentityWithRemoteAddress( const void* address ) const
    {
        return xEntity( address, getTag() );
    };

    void xMeshQueryInterface::printMesh() const
    {
        for ( xVertexIterator it = beginVertex(); it != endVertex(); ++it ) { ( *it ).print(); }
        for ( xEdgeIterator it = beginEdge(); it != endEdge(); ++it ) { ( *it ).print(); }
        for ( xFaceIterator it = beginFace(); it != endFace(); ++it ) { ( *it ).print(); }
        for ( xSolidIterator it = beginSolid(); it != endSolid(); ++it ) { ( *it ).print(); }
    };

    xEntityIterator xMeshQueryInterface::begin( int what ) const
    {
        switch ( what ) {
        case 0:
            return xEntityIterator( beginVertex() );
        case 1:
            return xEntityIterator( beginEdge() );
        case 2:
            return xEntityIterator( beginFace() );
        case 3:
            return xEntityIterator( beginSolid() );
        default:
            throw 6124;
        }
    };
    xEntityIterator xMeshQueryInterface::end( int what ) const
    {
        switch ( what ) {
        case 0:
            return xEntityIterator( endVertex() );
        case 1:
            return xEntityIterator( endEdge() );
        case 2:
            return xEntityIterator( endFace() );
        case 3:
            return xEntityIterator( endSolid() );
        default:
            throw 6125;
        }
    };

    void xMeshQueryInterface::exportGMSH( string filename ) const

    {

        string ext, text;
        auto   idx = filename.rfind( '.' );
        if ( idx > 0 ) {
            ext = filename.substr( idx );
            if ( ext == ".msh" ) {
                ofstream ofs( filename.c_str() );
                exportGMSH( ofs );
                ofs.close();
            } else {
                std::cout << "unknown extension " + ext + " in filename " << filename << " at line " << __LINE__ << " in file " << __FILE__ << std::endl;
                throw 6126;
            }
        } else {
            filename = filename + ".msh";
            ofstream ofs( filename.c_str() );
            exportGMSH( ofs );
            ofs.close();
        }
    };

    void xMeshQueryInterface::exportGMSH( ofstream& ofs ) const
    {
        ofs << "$NOE\n";
        ofs << getMeshSize( 0 ) << "\n";
        for ( auto it = beginVertex(); it != endVertex(); ++it ) {
            xVertex v = *it;

            ofs << v.getId() << " " << v.point()( 0 ) << " " << v.point()( 1 ) << " " << v.point()( 2 ) << "\n";
        }
        ofs << "$ENDNOE\n";
        int k = 0;
        int lev = 3;
        if ( ! getMeshSize( 3 ) ) lev = 2;
        if ( ! getMeshSize( 2 ) && ! getMeshSize( 3 ) ) lev = 1;
        for ( auto it = beginEdge(); it != endEdge(); ++it ) {
            xEntity edge = *it;
            if ( lev == 1 || edge.getClassification().dim() == 1 || edge.getClassification().tag() < 0 ) k++;
        }
        for ( auto it = beginFace(); it != endFace(); ++it ) {
            xEntity face = *it;
            if ( lev == 2 || face.getClassification().dim() == 2 || face.getClassification().tag() < 0 ) k++;
        }
        ofs << "$ELM\n";
        ofs << k + getMeshSize( 3 ) << "\n";
        k = 1;
        for ( auto it = beginEdge(); it != endEdge(); ++it ) {
            xEntity edge = *it;
            if ( lev == 1 || edge.getClassification().dim() == 1 )
                exportGmshEdge( ofs, edge, k, 0 );
        }
        for ( auto it = beginFace(); it != endFace(); ++it ) {
            xEntity face = *it;
            if ( lev == 2 || face.getClassification().dim() == 2 )
                exportGmshFace( ofs, face, k, 0 );
        }
        for ( auto it = beginSolid(); it != endSolid(); ++it )
            exportGmshSolid( ofs, ( *it ), k, 0 );

        ofs << "$ENDELM\n";
        ofs << "$FIN\n";
    }

    void exportGmshEdge( std::ostream& fs, const xEdge& e, int& k, int recur )
    {
        xVertex vx[2];
        for ( int i = 0; i < e.size( 0 ); i++ ) vx[i] = e.getVertex( i );
        for ( int i = 0; i < recur; i++ ) fs << " ";
        int nbrec = ( e.isAdjacencyCreated( 1 ) ) ? e.size( 1 ) : 0;
        fs << k++ << " 1 " << e.getClassification().tag() << " " << -nbrec << " "
           << 2 << " " << vx[0].getId() << " " << vx[1].getId() << "\n";
        for ( int i = 0; i < nbrec; i++ ) exportGmshEdge( fs, e.getEdge( i ), k, recur + 1 );
    }

    void exportGmshFace( std::ostream& fs, const xFace& e, int& k, int recur )
    {
        xVertex vx[4];
        for ( int i = 0; i < e.size( 0 ); i++ ) vx[i] = e.getVertex( i );
        for ( int i = 0; i < recur; i++ ) fs << " ";
        int nbrec = ( e.isAdjacencyCreated( 2 ) ) ? e.size( 2 ) : 0;
        int typ = ( e.size( 0 ) == 3 ) ? 2 : 3;
        if ( typ == 2 )
            fs << k++ << " " << typ << " " << e.getClassification().tag() << " "
               << e.getClassification().tag() << " "
               << 3 << " " << vx[0].getId() << " " << vx[1].getId() << " "
               << vx[2].getId() << "\n";
        else
            fs << k++ << " " << typ << " " << e.getClassification().tag()
               << " " << -nbrec << " "
               << 4 << " " << vx[0].getId() << " " << vx[1].getId()
               << " " << vx[2].getId() << " " << vx[3].getId() << "\n";
        for ( int i = 0; i < nbrec; i++ ) exportGmshFace( fs, e.getFace( i ), k, recur + 1 );
    }

    void exportGmshSolid( std::ostream& fs, const xSolid& e, int& k, int recur )
    {
        xVertex vx[8];
        for ( int i = 0; i < e.size( 0 ); i++ ) vx[i] = e.getVertex( i );
        for ( int i = 0; i < recur; i++ ) fs << " ";
        int nbrec = ( e.isAdjacencyCreated( 3 ) ) ? e.size( 3 ) : 0;
        switch ( e.getType() ) {
        case eType::PRISM:
            fs << k++ << " " << 6 << " " << e.getClassification().tag()
               << " " << -nbrec << " " << 6 << " " << vx[0].getId() << " " << vx[1].getId() << " " << vx[2].getId() << " " << vx[3].getId() << " " << vx[4].getId() << " " << vx[5].getId() << "\n";
            break;
        case eType::HEX:
            fs << k++ << " " << 5 << " " << e.getClassification().tag()
               << " " << -nbrec << " " << 8 << " " << vx[0].getId() << " " << vx[1].getId() << " " << vx[2].getId() << " " << vx[3].getId() << " " << vx[4].getId() << " " << vx[5].getId() << " " << vx[6].getId() << " " << vx[7].getId() << "\n";
            break;
        case eType::TET:
            fs << k++ << " " << 4 << " " << e.getClassification().tag()
               << " " << e.getClassification().tag() << " " << 4 << " " << vx[0].getId() << " " << vx[1].getId() << " " << vx[2].getId() << " " << vx[3].getId() << "\n";
            break;
        default: {
            std::cout << "Error in file " << __FILE__ << ":" << __LINE__
                      << " entity of type " << e.getType() << " not handled in switch."
                      << std::endl;
            throw 6129;
            break;
        }
        }
        for ( int i = 0; i < nbrec; i++ ) exportGmshSolid( fs, e.getSolid( i ), k, recur + 1 );
    }

} // namepsace xmeshinterface
} // namespace xinterface
