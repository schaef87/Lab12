//--------------------------------------------------------------------
//
//  Laboratory 12                                    WeightedGraph.cpp
//
//
//--------------------------------------------------------------------

#ifndef WEIGHTEDGRAPH_CPP
#define WEIGHTEDGRAPH_CPP

using namespace std;

#include "WeightedGraph.h"
#include "config.h"


//--------------------------------------------------------------------

WeightedGraph:: WeightedGraph ( int maxNumber )

// Creates an empty graph. Allocates enough memory for maxNumber
// vertices (defaults to defMaxGraphSize).

  : maxSize(maxNumber),
    size(0)
{
    ___________ = new Vertex [ maxSize ];  // <------------------------- fill in the member variable name

    ___________ = new int [ maxSize*maxSize ];  // <--------------------- fill in the member variable name
	                                                                      // NOTE: using 1 dimensional array to store the 2 dimensional
                                                                          // matrix (as describe on bottom of pg. 150)
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

WeightedGraph:: WeightedGraph ( const WeightedGraph& other )

// Creates an empty graph. Allocates enough memory for maxNumber
// vertices (defaults to defMaxGraphSize).

  : maxSize(other.maxSize),
    size(other.size)
{
    vertexList = new Vertex [ maxSize ];
    for( int i=0; i<size; i++ ) {
	vertexList[i] = other.vertexList[i];
    }

    adjMatrix = new int [ maxSize*maxSize ];
    for( int row=0; row<size; row++ ) {
	for( int col=0; col<size; col++ ) {
            adjMatrix[maxSize*row+col] = other.adjMatrix[maxSize*row+col];
	}
    }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

WeightedGraph& WeightedGraph:: operator= ( const WeightedGraph& other )

// Copies from another WeightedGraph.

{
    if( this == &other ) {	// Self-assignment attempt!
        return *this;
    }

    // If arrays are not large enough, delete current and create new
    if( maxSize < other.maxSize ) {
	delete[] vertexList;
	delete[] adjMatrix;

        vertexList = new Vertex [ other.maxSize ];
        adjMatrix = new int [ other.maxSize*other.maxSize ];
    }

    maxSize = other.maxSize;
    size = other.size;

    for( int i=0; i<size; i++ ) {
	vertexList[i] = other.vertexList[i];
    }

    for( int row=0; row<size; row++ ) {
	for( int col=0; col<size; col++ ) {
            adjMatrix[maxSize*row+col] = other.adjMatrix[maxSize*row+col];
	}
    }
}

//--------------------------------------------------------------------

WeightedGraph:: ~WeightedGraph ()

// Frees the memory used by a graph.

{
    delete [] vertexList;
    delete [] adjMatrix;
}

//--------------------------------------------------------------------

void WeightedGraph:: insertVertex ( const Vertex& newVertex ) throw ( logic_error )

// Inserts newVertex into a graph. If a vertex with the same label
// as newVertex already exists in the graph, then updates that
// vertex's data with newVertex's data.

{
    int pos,   // Position of the vertex in the vertex list
        j;     // Loop counter

    // Requires that graph is not full
    if ( size == maxSize )
        throw logic_error("graph is full");

    pos = getIndex(newVertex.getLabel());
    vertexList[ ___ ] = newVertex;  // <------------------------- fill in the array index

// debug
//cout << "newVertex.getLabel = " << newVertex.getLabel()
//     << ", pos = " << pos << endl;

    if ( pos == size )            // New vertex
    {
       size++;
       for ( j = 0 ; j < size ; j++ )
       {
           setEdge(size-1, j, INFINITE_EDGE_WT);
           setEdge(j, size-1, INFINITE_EDGE_WT);
       }
    }
}

//--------------------------------------------------------------------

void WeightedGraph:: insertEdge ( const string& v1, const string& v2, int wt )
    throw ( logic_error )

// Insert an edge with the specified weight between vertices
// v1 and v2.

{
    int idx_v1 = getIndex(v1),   // Index of vertex v1
        idx_v2 = getIndex(v2);   // Index of vertex v2

    //cout << "Size = " << size << ", idx_v1 = " << idx_v1 << ", idx_v2 = " << idx_v2 << endl;
    if ( idx_v1 >= size )
        throw logic_error(string("idx_v1 is out of range. v1=")+v1);

    if ( idx_v2 >= size )
        //throw logic_error("idx_v2 is out of range");
        throw logic_error(string("idx_v2 is out of range. v2=")+v2);

    // set edge weight in adj. matrix using helper function
	_______(idx_v1,idx_v2,wt);  // <------------------------------------ fill in the helper function name
    _______(idx_v2,idx_v1,wt);  // <------------------------------------ fill in the helper function name
}

//--------------------------------------------------------------------

bool WeightedGraph:: retrieveVertex ( const string& v, Vertex &vData ) const

// Searches a graph for vertex v. If the vertex is found, then copies
// the vertex data to vData and returns 1. Otherwise, returns 0 with
// vData undefined.

{
    int pos;      // Position of the vertex in the vertex list
    bool result;  // Result returned

    pos = getIndex(v);
    if ( pos != size )
    {
       vData = vertexList[pos];
       result = true;
    }
    else
       result = false;

    return result;
}

//--------------------------------------------------------------------

bool WeightedGraph:: getEdgeWeight ( const string& v1, const string& v2, int &wt ) const

    throw ( logic_error )

// If there is an edge connecting vertices v1 and v2, then returns 1
// with wt returning the weight of the edge. Otherwise, returns 0
// with wt undefined.

{
    int idx_v1 = ____________ ;   // Index of vertex v1    // <------------------------- complete this line
    int idx_v2 = ____________ ;   // Index of vertex v2    // <------------------------- complete this line
	                                                                                 // NOTE: this function works similarly
																					 // to the insertEdge function above

    if ( idx_v1 >= size )
        throw logic_error("idx_v1 is out of range");

    if ( idx_v2 >= size )
        throw logic_error("idx_v2 is out of range");

    wt = ______________  ;  // Get weight of edge using helper function  <------------------------- complete this line
    return ( wt != INFINITE_EDGE_WT );
}

//--------------------------------------------------------------------

void WeightedGraph:: removeVertex ( const string& v ) throw ( logic_error )

// Removes vertex v from a graph.

{
    int idx_v = getIndex(v),   // Index of vertex v
        j, k;                  // Loop counters

    if ( idx_v >= size )                     // Vertex is required to exist!
        throw logic_error("idx_v is out of range");

    for ( j = idx_v ; j < size-1 ; j++ )     // Shift vertex list up
        vertexList[j] = vertexList[j+1];

    for ( j = idx_v ; j < size-1 ; j++ )     // Shift rows up
        for ( k = 0 ; k < size ; k++ )
            setEdge(j,k, getEdge(j+1,k));

    for ( j = idx_v ; j < size-1 ; j++ )     // Shift columns left
        for ( k = 0 ; k < size ; k++ )
            setEdge(k,j, getEdge(k,j+1));

    size--;
}

//--------------------------------------------------------------------

void WeightedGraph:: removeEdge ( const string& v1, const string& v2 ) throw ( logic_error )

// Removes the edge between vertices v1 and v2 from a graph.

{
    int idx_v1 = getIndex(v1),   // Index of vertex v1
        idx_v2 = getIndex(v2);   // Index of vertex v2

    if ( idx_v1 >= size )
        throw logic_error("idx_v1 is out of range");

    if ( idx_v2 >= size )
        throw logic_error("idx_v2 is out of range");

    setEdge(idx_v1,idx_v2, INFINITE_EDGE_WT);
    setEdge(idx_v2,idx_v1, INFINITE_EDGE_WT);
}

//--------------------------------------------------------------------

void WeightedGraph:: clear ()

// Removes all the vertices and edges from a graph.

{
    size = 0;
}

//--------------------------------------------------------------------

bool WeightedGraph:: isEmpty () const

// Returns 1 if a graph is empty. Otherwise, returns 0.

{
    return ( size == 0 );
}

//--------------------------------------------------------------------

bool WeightedGraph:: isFull () const

// Returns 1 if a graph is full. Otherwise, returns 0.

{
    return ( size == maxSize );
}

//--------------------------------------------------------------------

void WeightedGraph:: showStructure () const

// Outputs a graph's vertex list and adjacency matrix. This operation
// is intended for testing/debugging purposes only.

{
    if ( size == 0 ) {
       cout << "Empty graph" << endl;
    } else {
       cout << endl << "Vertex list : " << endl;
       for ( int row = 0 ; row < size ; row++ ) {
           cout << row << '\t' << vertexList[row].getLabel();
#if LAB12_TEST2
	   cout << '\t' << vertexList[row].getColor();
#endif
	   cout << endl;
       }

       cout << endl << "Edge matrix : " << endl << '\t';
       for ( int col = 0 ; col < size ; col++ )
           cout << col << '\t';
       cout << endl;
       for ( int row = 0 ; row < size ; row++ )
       {
           cout << row << '\t';
           for ( int col = 0 ; col < size ; col++ )
           {
               int wt = getEdge(row,col);
               if ( wt == INFINITE_EDGE_WT )
                  cout << "- \t";
               else
                  cout << wt << '\t';
           }
           cout << endl;
       }
    }
}

#if LAB12_TEST1		// computePaths
class pathMatrix {
  public:
    pathMatrix(int initSize) : size(initSize) {
	matrix = new int[initSize * initSize];
    }

    ~pathMatrix() {
	delete matrix;
    }

    int getPath(int row, int col) const {
	return matrix[size * row + col];
    }

    void setPath(int row, int col, int wt) {
	matrix[size * row + col] = wt;
    }

  private:
    int *matrix;
    int size;
};

void WeightedGraph::showShortestPaths() const {
    pathMatrix pm(size);

    for ( int j = 0 ; j < size ; j++ )           // Copy edge matrix
        for ( int k = 0 ; k < size ; k++ )
            pm.setPath(j,k, getEdge(j,k));

    for ( int j = 0 ; j < size ; j++ )           // Set diagonal to 0
        pm.setPath(j,j, 0);

    for ( int m = 0 ; m < size ; m++ )           // Compute paths
        for ( int j = 0 ; j < size ; j++ )
            for ( int k = 0 ; k < size ; k++ )
                if ( long(pm.getPath(j,m)) + pm.getPath(m,k) < pm.getPath(j,k) )
                   //pm.setPath(j,k, (_______________________________)))   // <------------ After your WeightedGraph is working, uncomment this line and fill in
			         ;                                                                   // the blank using Floyd's algorithm described in the middle of page. 157

    cout << endl << "Path matrix : " << endl << '\t';
    for ( int col = 0 ; col < size ; col++ )
       cout << col << '\t';
    cout << endl;
    for ( int row = 0 ; row < size ; row++ )
    {
       cout << row << '\t';
       for ( int col = 0 ; col < size ; col++ )
       {
	   int wt = pm.getPath(row,col);
	   if ( wt == INFINITE_EDGE_WT )
	      cout << "- \t";
	   else
	      cout << wt << '\t';
       }
       cout << endl;
    }
}
#endif



//--------------------------------------------------------------------
//
//  Facilitator functions
//

int WeightedGraph:: getIndex ( const string& v ) const

// Returns the adjacency matrix index for vertex v. Returns size if
// the vertex does not exist.

{
    int j;  // Loop counter

    for ( j = 0 ;
          j < size  &&  (vertexList[j].getLabel() != v);
          j++ );
// debug
//cout << "getIndex: vertex = " << v << ", pos = " << j << endl;

    return j;
}

//--------------------------------------------------------------------

int WeightedGraph:: getEdge ( int row, int col ) const

// Gets adjMatrix[row][col].

{
    return adjMatrix[ ______________ ];  // <---------------------------------------------- fill in the array index
}

void WeightedGraph:: setEdge ( int row, int col, int wt )

// Gets adjMatrix[row][col].

{
    adjMatrix[maxSize*row+col] = wt;
}

#endif	// LAB12_TEST3

#endif	// #ifndef WEIGHTEDGRAPH_CPP
