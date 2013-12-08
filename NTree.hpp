#ifndef NTREE_H
#define NTREE_H

#include "General.hpp"
#include "Vec3.hpp"
#include "Box.hpp"

// A class templated to create
// DIM = 1: BinaryTrees    Vec3's are (x,0,0)
// DIM = 2: QuadTrees      Vec3's are (x,y,0)
// DIM = 3: OctTrees       Vec3's are (x,y,z)


template <int DIM>
class NTree
{
  double rootSize;             // The size of the root box
  Vec3 pmin;                   // Handle to the bottom left of the root box

  int nLevels;                 // The number of levels: 0 = root, nLevel = leaf
  int nBoxes;

  // Should sparsify this if very few boxes are actually required
  vector<Box*> boxes;          // The boxes in this tree

 public:

  const static int BRANCH = 1 << DIM;   // The branching factor = 2^DIM

  // Constructor of an NTree with levels [0,L]
  NTree( vector<Vec3>& p, int L )
      : nLevels(L),
        nBoxes( ((1 << (DIM*(L+1)))-1)/(BRANCH-1) ),
        boxes( nBoxes )
  {
    int N = p.size();

    // Find normalization to bounding box [0,1]^DIM
    pmin = Vec3( INFINITY, INFINITY, INFINITY );
    Vec3 pmax = -pmin;
    for( int k = 0; k < N; ++k ) {
      Vec3& t = p[k];

      if( DIM <= 1 ) p[k].y = 0;
      if( DIM <= 2 ) p[k].z = 0;

      if( t.x < pmin.x ) pmin.x = t.x;
      if( t.x > pmax.x ) pmax.x = t.x;
      if( t.y < pmin.y ) pmin.y = t.y;
      if( t.y > pmax.y ) pmax.y = t.y;
      if( t.z < pmin.z ) pmin.z = t.z;
      if( t.z > pmax.z ) pmax.z = t.z;
    }

    if( DIM <= 2 ) assert(pmin.z == pmax.z);
    if( DIM <= 1 ) assert(pmin.y == pmax.y);

    // The sidelength of the largest (encompassing) box
    rootSize = (1+1e-14)*max(max(pmax.x-pmin.x,pmax.y-pmin.y),pmax.z-pmin.z);

    // Remember: Using ints as bit array only allows finite L
    int max_bits = 30;
    int coord_bit_count = max_bits/DIM;
    assert( nLevels <= coord_bit_count );

    // Define and order points in the normalized box
    double scale = (1 << coord_bit_count)/rootSize;
    vector< pair<int,int> > mPair(N);
    for( int k = 0; k < N; ++k ) {
      Vec3& t = p[k];

      int pNum = interleave( (int)(scale*(t.x-pmin.x)),
                             (int)(scale*(t.y-pmin.y)),
                             (int)(scale*(t.z-pmin.z)) );

      mPair[k] = pair<int,int>( pNum, k );
    }
    sort( mPair.begin(), mPair.end() );

    // Create lowest level of boxes
    Box* currentBox = NULL;
    int currentBoxN = -1;
    // For each (ordered) point
    for( int k = 0; k < N; ++k ) {
      pair<int,int>& currentPoint = mPair[k];

      // Determine the box number it belongs in
      int boxN = parent( currentPoint.first, coord_bit_count - nLevels);
      int pointN = currentPoint.second;

      // If this is a new box number, create and add the box to the level
      if( boxN != currentBoxN ) {
        currentBox = new Box( boxN, boxCenter(boxN, nLevels) );
        boxes[ getBoxIndex(boxN, nLevels) ] = currentBox;
        currentBoxN = boxN;
      }

      // Add the point index to the box it is contained in
      currentBox->addPointIndex( pointN );
    }

    // Create the boxes for the rest of the levels
    for( int L = nLevels; L > 0; --L ) {
      Box* currentParent = NULL;
      int currentParentN = -1;

      // For each box
      for( int n = 0; n < (1 << (DIM*L)); ++n ) {
        Box* b = boxes[ getBoxIndex(n,L) ];
        if( b == NULL ) continue;

        // Determine the parent box number it belongs in
        int parentN = parent( b->n );

        // If this is a new box number, create and add the box to the level
        if( parentN != currentParentN ) {
          currentParent = new Box( parentN, boxCenter( parentN, L-1 ) );
          boxes[ getBoxIndex(parentN, L-1) ] = currentParent;
          currentParentN = parentN;
        }

        // Set the parent pointer
        b->parent = currentParent;
      }
    }
  }

  // Destructor
  ~NTree() {
    // Delete Boxes
    for( int n = 0; n < nBoxes; ++n ) {
      delete boxes[n];
    }
  }

  inline int numLevels()
  {
    return nLevels;
  }

  inline int getBoxIndex( int n, int L )
  {
    assert( L <= nLevels );
    assert( n < (1 << (DIM*L)) );
    //cerr << n << ", " << L << ": " << ((1 << (DIM*L))-1)/(BRANCH-1) + n << " -- " << nBoxes << endl;
    return ((1 << (DIM*L))-1)/(BRANCH-1) + n;
  }

  inline Box* getBox( int n, int L )
  {
    return boxes[ getBoxIndex(n,L) ];
  }

  // Note that L = 0 is the root
  inline double getBoxSize( int L )
  {
    return rootSize / (1 << L);
  }

  inline Vec3 boxCenter( int n, int L )
  {
    Vec3 p = deinterleave(n);

    if( DIM >= 1 ) p.x += 0.5;
    if( DIM >= 2 ) p.y += 0.5;
    if( DIM >= 3 ) p.z += 0.5;

    p *= getBoxSize(L);
    p += pmin;
    return p;
  }

  // The kth-order parent of box n
  inline int parent( int n, int k = 1 )
  {
    return n >> (DIM*k);
  }

  // The mth kth-order child of box n
  inline int child( int n, int m = 0, int k = 1 )
  {
    assert( m < (1 << (DIM*k)) );
    return m + (n << (DIM*k));
  }

  inline void addPoint( const Vec3& p ) {}

  inline void removePoint( int index ) {}

  friend ostream& operator<<(ostream& os, NTree& t)
  {
    int nLevels = t.numLevels();
    for( int n = 0; n < (1 << (DIM*nLevels)); ++n ) {
      Box* b = t.getBox(n,nLevels);
      if( b == NULL ) continue;

      for( int l = nLevels; l >= 0; --l ) {
        os << t.parent(b->n,l) << "\t";
      }
      os << b->center << endl;
    }
    return os;
  }

 private:

  inline int interleave( int x, int y, int z )
  {
    int n = 0;

    for( int i = -1; x != 0 || y != 0 || z != 0; ) {
      if( DIM >= 3 ) { n |= (z & 1) << ++i;  z >>= 1; }
      if( DIM >= 2 ) { n |= (y & 1) << ++i;  y >>= 1; }
      if( DIM >= 1 ) { n |= (x & 1) << ++i;  x >>= 1; }
    }

    return n;
  }

  inline Vec3 deinterleave( int n )
  {
    int x = 0, y = 0, z = 0;

    for( int i = 0; n != 0; ++i ) {
      if( DIM >= 3 ) { z |= (n & 1) << i;  n >>= 1; }
      if( DIM >= 2 ) { y |= (n & 1) << i;  n >>= 1; }
      if( DIM >= 1 ) { x |= (n & 1) << i;  n >>= 1; }
    }

    return Vec3(x,y,z);
  }

  inline vector<int> siblings( int n )
  {
    vector<int> siblings(BRANCH);
    siblings[0] = child(parent(n));
    for( int k = 1; k < BRANCH; ++k ) {
      siblings[k] = siblings[k-1] + 1;
    }
    return siblings;
  }

 private:
  // Disable Copy and Assignment
  NTree( const NTree& S ) {}
  void operator=( const NTree& S ) {}
};


#endif
