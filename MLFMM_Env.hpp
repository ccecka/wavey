#ifndef MLFMM_Env_CPP
#define MLFMM_Env_CPP

#include "NTree.hpp"
#include "Level.hpp"

/*
 * A class to store MLFMM data which should be reused between iterations
 */

#define STAGE_TIMES 0


template <typename FUNC, int DIM>
class MLFMM_Env
{
 public:
  FUNC K;

  vector<Vec3> sourcefield;

  NTree<DIM> tree;
  vector<Level<DIM>> levels;

  // Constructor
  MLFMM_Env( FUNC K_, vector<Vec3>& p,
             int levels_ = 4, double eps = 1e-4 )
      : K(K_),
        sourcefield( p ),       // Copy the charge positions
        tree( p, levels_ ),     // Create the tree
        levels( levels_+1 )     // Initialize vector of levels
  {
    cerr << tree << endl;

    int nLevels = tree.numLevels();

    // Construct each Level and give it the boxes it owns
    for( int L = 2; L <= nLevels; ++L ) {
      levels[L] = Level<DIM>( tree.getBoxSize(L) );
      for( int n = 0; n < (1 << (DIM*L)); ++n ) {
        Box* b = getBox(n,L);
        if( b != NULL )
          getLevel(L).addBox(b);
      }
    }

    // Recursively build transfer and close lists from rootbox
    defineTransferAndClose(getBox(0,0), 0);

    cerr << "Transfer and Close Lists Built" << endl;

    // Compute the quadrature for each level of the tree
    for( int L = 2; L <= nLevels; ++L ) {
      getLevel(L).defineQuadrature( K, eps );
    }

    cout << "Quadrature Defined" << endl;

    // Compute the transfer functions for each level of the tree
    for( int L = 2; L <= nLevels; ++L ) {
      getLevel(L).defineTransfers();
    }

    cout << "Transfer Functions Defined" << endl;

    // Compute the translation operators for each level of the tree
    for( int L = 2; L < nLevels; ++L ) {
      getLevel(L).defineTranslations();
    }

    cout << "Translation Functions Defined" << endl;

    // Define the interpolators and anterpolators required for each level
    for( int L = 2; L <= nLevels; ++L ) {
      if( L != 2 )
        getLevel(L).defineInterp( getLevel(L-1).getQuad() );
      if( L != nLevels )
        getLevel(L).defineAnterp( getLevel(L+1).getQuad() );
    }

    cout << "Reterpolators Defined" << endl;

    cout << "Initialization Done" << endl;
  }

  // Destructor
  ~MLFMM_Env() {}

  // Computes   omega = A*psi
  void execute( const vector<complex>& psi, vector<complex>& omega )
  {
#if STAGE_TIMES
    StopWatch timer;
    timer.start();
#endif

    int nLevels = tree.numLevels();

    // Compute lowest level multipole
    Level<DIM>& leafLevel = getLevel(nLevels);
    leafLevel.zeroFields();
    for( BoxIter bi = leafLevel.boxbegin(); bi != leafLevel.boxend(); ++bi ) {
      Box* b = *bi;

      const Vec3& center = b->center;

      const Box::pointIter piEnd = b->pointIndex.end();
      for( Box::pointIter pi = b->pointIndex.begin(); pi != piEnd; ++pi ) {
        int index = *pi;
        const Vec3 r = center - sourcefield[index];
        cout << "Accumulating " << index << " into Box " << b->n << endl;
        Translation_Function::add( r, psi[index], b->getMultipole() );
      }
    }

#if STAGE_TIMES
    cout << "Lowest Level Multipoles Done " << timer.stop() << endl;
#endif

    // Upward Pass - Compute the multipoles
    for( int L = nLevels; L >= 3; --L ) {

      Level<DIM>& level = getLevel(L);
      Level<DIM>& pLevel = getLevel(L-1);
      pLevel.zeroFields();
      NFunction& Ms = pLevel.getScratch();

      // For all boxes at level L
      for( BoxIter bi = level.boxbegin(); bi != level.boxend(); ++bi ) {
        Box* b = *bi;

        // Interpolate to the parent level L-1
        level.getInterp().apply( b->getMultipole(), Ms );

        // Multiply by the transfer function and accumulate into parent
        b->parent->getMultipole().addProduct( Ms, pLevel.getTranslationUp(b) );
      }
    }

#if STAGE_TIMES
    cout << "Upward Pass Done " << timer.stop() << endl;
#endif

    // Do the Transfers
    for( int L = nLevels; L >= 2; --L ) {
      getLevel(L).applyTransferFunctions();
    }

#if STAGE_TIMES
    cout << "Transfer List Applied " << timer.stop() << endl;
#endif

    for( int L = 3; L <= nLevels; ++L ) {

      Level<DIM>& level = getLevel(L);
      Level<DIM>& pLevel = getLevel(L-1);
      NFunction& Lb = level.getScratch();
      NFunction& LB = pLevel.getScratch();

      // For all boxes at level L
      for( BoxIter bi = level.boxbegin(); bi != level.boxend(); ++bi ) {
        Box* b = *bi;

        // Multiply parent's local expansion by translation function
        LB.setProduct( b->parent->getLocal(), pLevel.getTranslationDown(b) );

        // Anterpolate
        pLevel.getAnterp().apply( LB, Lb );

        // Accumulate into the box's local expansion
        b->getLocal().add( Lb );
      }
    }

#if STAGE_TIMES
    cout << "Downward Pass Done " << timer.stop() << endl;
#endif

    // Zero omega
    omega.assign( omega.size(), 0 );

    // Integrate and add close interactions

    // For all the boxes at the lowest level, integrate local
    NFunction& LE = leafLevel.getScratch();           // Scratch space
    for( BoxIter bi = leafLevel.boxbegin(); bi != leafLevel.boxend(); ++bi ) {
      Box* b = *bi;
      Vec3& center = b->center;
      // for all interior points
      const Box::pointIter nEnd = b->pointIndex.end();
      for( Box::pointIter ni = b->pointIndex.begin(); ni != nEnd; ++ni ) {
        int n = *ni;
        const Vec3& pn = sourcefield[n];
        Vec3 r = pn - center;
        //cout << "Integrating for " << n << " from box " << b->n << endl;
        Translation_Function::times( b->getLocal(), r, LE );
        omega[n] += LE.integrate();

        // Add contributions from other points inside this box
        for( Box::pointIter mi = b->pointIndex.begin(); mi != ni; ++mi ) {
          int m = *mi;
          double r = pn.dist(sourcefield[m]);
          complex E = K(r);
          omega[m] += psi[n] * E;
          omega[n] += psi[m] * E;
        }
      }
    }

#if STAGE_TIMES
    cout << "Local Integration Done " << timer.stop() << endl;
#endif

    leafLevel.applyClose(K, psi, sourcefield, omega);

#if STAGE_TIMES
    cout << "Close List Applied " << timer.stop() << endl;
#endif
  }

  // Should be private...
  inline Level<DIM>& getLevel( int L ) { return levels[L]; }
  inline Box* getBox( int n, int L ) { return tree.getBox(n,L); }

 private:

  // Main algorithm to determine transfer pairs and close pairs
  // Two boxes are transfers if their parents are neighbors and they aren't
  // Two boxes are close if they are leaf level and are neighbors
  inline void defineTransferAndClose( Box* b, int L )
  {
    // For all the children of b
    for( int k1 = 0; k1 < tree.BRANCH; ++k1 ) {
      // Get the k1th child of b
      Box* child1 = getBox( tree.child(b->n, k1), L+1 );
      if( child1 == NULL ) continue;

      // For all the other children of b
      for( int k2 = k1+1; k2 < tree.BRANCH; ++k2 ) {
        // Get the k2th child of b
        Box* child2 = getBox( tree.child(b->n, k2), L+1 );
        if( child2 == NULL ) continue;

        // Call the neighbor routine on each pair of children
        defineTransferAndClose( child1, child2, L+1 );
      }

      // If the child is not at the leaf level, recurse on the child
      if( L+1 != tree.numLevels() )
        defineTransferAndClose( child1, L+1 );
    }
  }

  // Secondary algorithm to determine transfer pairs and close pairs
  inline void defineTransferAndClose( Box* b1, Box* b2, int L )
  {
    // If leaf level, these boxes are a close pair
    if( L == tree.numLevels() ) {
      getLevel(L).addClose( b1, b2 );
      return;
    }

    // For all the children of b1
    for( int k1 = 0; k1 < tree.BRANCH; ++k1 ) {
      // Get the k1th child of b1
      Box* child1 = getBox( tree.child(b1->n, k1), L+1 );
      if( child1 == NULL ) continue;

      // For all the children of b2
      for( int k2 = 0; k2 < tree.BRANCH; ++k2 ) {
        // Get the k2th child of b2
        Box* child2 = getBox( tree.child(b2->n, k2), L+1 );
        if( child2 == NULL ) continue;

        // Determine the transfer vector between the two children
        Vec3 r0 = child1->center - child2->center;
        if( r0.mag() > 1.8 * tree.getBoxSize(L+1) ) {
          // These two are not neighbors so they are a transfer pair
          getLevel(L+1).addTransfer( child1, child2 );
        } else {
          // These two are neighbors so recurse
          defineTransferAndClose( child1, child2, L+1 );
        }
      }
    }
  }

};

#endif
