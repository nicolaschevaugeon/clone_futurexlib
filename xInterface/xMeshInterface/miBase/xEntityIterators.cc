#include "xEntityIterators.h"
#include "xMeshQueryInterface.h"
#include <iostream>

namespace xinterface {

namespace xmeshinterface {

    using namespace std;

    // --------VERTEX -------

    xVertex xVertexIterator::operator*() const
    {
        return pQuery->fromVertexIteratorIdentifier( iterator_identifier );
    }
    xVertexIterator& xVertexIterator::operator++()
    {
        iterator_identifier = pQuery->nextVertexIteratorIdentifier( iterator_identifier );
        return *this;
    }
    xVertexIterator xVertexIterator::operator++( int )
    {
        xVertexIterator old( *this );
        ++( *this );
        return old;
    }
    bool xVertexIterator::operator!=( const xVertexIterator& in ) const
    {
        if ( pQuery != in.pQuery ) return true;
        return ! pQuery->equalVertexIteratorIdentifier( iterator_identifier, in.iterator_identifier );
    }
    bool xVertexIterator::operator==( const xVertexIterator& in ) const
    {
        if ( pQuery != in.pQuery ) return false;
        return pQuery->equalVertexIteratorIdentifier( iterator_identifier, in.iterator_identifier );
    }

    // --------EDGE -------

    xEdge xEdgeIterator::operator*() const
    {
        return pQuery->fromEdgeIteratorIdentifier( iterator_identifier );
    }
    xEdgeIterator& xEdgeIterator::operator++()
	{
        iterator_identifier = pQuery->nextEdgeIteratorIdentifier( iterator_identifier );
        return *this;
	}
    xEdgeIterator xEdgeIterator::operator++( int )
	{
        xEdgeIterator old( *this );
        ++( *this );
        return old;
	}
    bool xEdgeIterator::operator!=( const xEdgeIterator& in ) const
    {
        if ( pQuery != in.pQuery ) return true;
        return ! pQuery->equalEdgeIteratorIdentifier( iterator_identifier, in.iterator_identifier );
	}
    bool xEdgeIterator::operator==( const xEdgeIterator& in ) const
	{
        if ( pQuery != in.pQuery ) return false;
        return pQuery->equalEdgeIteratorIdentifier( iterator_identifier, in.iterator_identifier );
	}

    // --------FACE -------
    xFace xFaceIterator::operator*() const
    {
        return pQuery->fromFaceIteratorIdentifier( iterator_identifier );
    }
    xFaceIterator& xFaceIterator::operator++()
    {
        iterator_identifier = pQuery->nextFaceIteratorIdentifier( iterator_identifier );
        return *this;
    }
    xFaceIterator xFaceIterator::operator++( int )
    {
        xFaceIterator old( *this );
        ++( *this );
        return old;
    }
    bool xFaceIterator::operator!=( const xFaceIterator& in ) const
    {
        if ( pQuery != in.pQuery ) return true;
        return ! pQuery->equalFaceIteratorIdentifier( iterator_identifier, in.iterator_identifier );
    }
    bool xFaceIterator::operator==( const xFaceIterator& in ) const
    {
        if ( pQuery != in.pQuery ) return false;
        return pQuery->equalFaceIteratorIdentifier( iterator_identifier, in.iterator_identifier );
    }

    // --------SOLID -------
    xSolid xSolidIterator::operator*() const
	{
        return pQuery->fromSolidIteratorIdentifier( iterator_identifier );
	}
    xSolidIterator& xSolidIterator::operator++()
	{
        iterator_identifier = pQuery->nextSolidIteratorIdentifier( iterator_identifier );
        return *this;
	}
    xSolidIterator xSolidIterator::operator++( int )
	{
        xSolidIterator old( *this );
        ++( *this );
        return old;
	}
    bool xSolidIterator::operator!=( const xSolidIterator& in ) const
	{
        if ( pQuery != in.pQuery ) return true;
        return ! pQuery->equalSolidIteratorIdentifier( iterator_identifier, in.iterator_identifier );
	}
    bool xSolidIterator::operator==( const xSolidIterator& in ) const
	{
        if ( pQuery != in.pQuery ) return false;
        return pQuery->equalSolidIteratorIdentifier( iterator_identifier, in.iterator_identifier );
    }

    // -------- ENTITY -------
    xEntity xEntityIterator::operator*() const
    {
        switch ( type ) {
        case VERTEXITERATOR: {
            xVertexIterator it( iter );
            return xEntity( *it );
            break;
        }
        case EDGEITERATOR: {
            xEdgeIterator it( iter );
            return xEntity( *it );
            break;
        }
        case FACEITERATOR: {
            xFaceIterator it( iter );
            return xEntity( *it );
            break;
        }
        case SOLIDITERATOR: {
            xSolidIterator it( iter );
            return xEntity( *it );
            break;
        }
        default:
            throw 49654;
            return xEntity();
        }
    };

    xEntityIterator& xEntityIterator::operator++()
	{

        switch ( type ) {
        case VERTEXITERATOR: {
            xVertexIterator it( this->iter );
            ++it;
            iter = it;
            return *this;
            break;
        }
        case EDGEITERATOR: {
            xEdgeIterator it( this->iter );
            ++it;
            iter = it;
            return *this;
            break;
        }
        case FACEITERATOR: {
            xFaceIterator it( this->iter );
            ++it;
            iter = it;
            return *this;
            break;
        }
        case SOLIDITERATOR: {
            xSolidIterator it( this->iter );
            ++it;
            iter = it;
            return *this;
            break;
        }
        default:
            throw 49654;
            return *this;
        }
    };

    xEntityIterator xEntityIterator::operator++( int )
	{
        xEntityIterator old( *this );
        ++( *this );
        return old;
    };

    bool xEntityIterator::operator!=( const xEntityIterator& in ) const
	{
        if ( type != in.type )
            return true;

        switch ( type ) {

        case VERTEXITERATOR: {
            xVertexIterator it( this->iter );
            xVertexIterator sh( in.iter );
            return it != sh;
            break;
        }
        case EDGEITERATOR: {
            xEdgeIterator it( this->iter );
            xEdgeIterator sh( in.iter );
            return it != sh;
            break;
        }
        case FACEITERATOR: {
            xFaceIterator it( iter );
            xFaceIterator sh( in.iter );
            return it != sh;
            break;
        }
        case SOLIDITERATOR: {
            xSolidIterator it( iter );
            xSolidIterator sh( in.iter );
            return it != sh;
            break;
        }
        default:
            throw 49654;
            return false;
        }
    };

    bool xEntityIterator::operator==( const xEntityIterator& in ) const
    {
        if ( type != in.type )
            return false;

        switch ( type ) {
        case VERTEXITERATOR: {
            xVertexIterator it( this->iter );
            xVertexIterator sh( in.iter );
            return it == sh;
            break;
        }
        case EDGEITERATOR: {
            xEdgeIterator it( this->iter );
            xEdgeIterator sh( in.iter );
            return it == sh;
            break;
        }
        case FACEITERATOR: {
            xFaceIterator it( this->iter );
            xFaceIterator sh( in.iter );
            return it == sh;
            break;
        }
        case SOLIDITERATOR: {
            xSolidIterator it( this->iter );
            xSolidIterator sh( in.iter );
            return it == sh;
            break;
        }
        default:
            throw 49654;
            return false;
        }
    };

} // namepsace xmeshinterface
} // namespace xinterface
